% Perform Snyder filtering for variable N(t) = N(0)exp(-xt) (reverse time)

% Assumptions and Modifications
% - added repeatability for M runs
% - exponential formulation for lam
% - used time rescaling from Palacios 2012 to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clear all

% Booleans to control functionality
M = 1;
plotindiv = 1;


% Parameters to store results
MSEpercratio = zeros(1, M);
param = zeros(1, M);
paramhat = zeros(1, M);
Cset = zeros(1, M);
rset = zeros(1, M);
MSE = zeros(1, M);

for ii = 1:M
    
    %% Generate the appropriate coalescent data for a defined N(t)
    
    % Set initial parameters: q0 is prior, N is true value, m is dimension
    % of parameter space, n is no. samples and m is total vector length
    m = 100;
    q0 = ones(1, m)/m;
    n = 1000;
    nData = n-1;
    
    % Set population parameters for exponential with N0 >> n
    xmax = 100;
    xmin = 0.001;
    xset = linspace(xmin, xmax, m);
    x = randsample(xset, 1);
    N0 = n*(10^3);
    
    % Get coalescent binomial factors for lineages
    nset = n:-1:2;
    fac = nset.*(nset-1)/2;
    
    % Simulate coalescent times for exponential using time rescaling
    tcoal = zeros(1, n);
    for i = 1:n-1
        % Get coalescent time for i+1 -> i lineages
        tcoal(i+1) = log(exp(x*tcoal(i)) + x*N0*exprnd(1)/fac(i))/x;
    end
    % Get waiting times
    twait = diff(tcoal);
    tmrca = tcoal(end);
    
    %% Filter the simulated event times from the appropriate Poisson process
    
    % Posterior vectors on events, qnoev is for non-updated event q's
    qev = zeros(nData+1, m);
    qnoev = qev;
    qev(1, :) = q0;
    tev = zeros(nData, 1);
    
    % Cell to save output of ODE solver and set options
    qset = cell(1, 1);
    tset = cell(1, 1);
    xhatset = cell(1, 1);
    xhatset2 = cell(1, 1);
    thetaset = cell(1, 1);
    lamset = cell(1, 1);
    Nset = cell(1, 1);
    options = odeset('NonNegative', 1:m);
    elemLen = zeros(1, nData);
    xev = zeros(1, nData);
    
    % Run Snyder filter across the time series
    for i = 1:nData
        % Obtain approproate binomial factor
        binfac = fac(i);
        
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnyExp1Param(ts, y, [], binfac, xset, N0),...
            [tcoal(i) tcoal(i+1)], qev(i, :)', options);
        
        % Normalise the posterior probabilities
        for j = 1:size(qsol, 1)
            qsol(j, :) = qsol(j, :)/(sum(qsol(j, :)));
        end
        
        % Assign the output values of time and posteriors
        qnoev(i, :) = qev(i, :);
        qset{i} = qsol;
        tset{i} = tsol;
        elemLen(i) = length(tsol);
        
        % Calculate lamt across time explicitly
        lamtdiag = zeros(length(tsol), m);
        for j = 1:length(tsol)
            % Calculate exponential inhomogeneous rate matrix diagonal
            lamtdiag(j, 1:m) = (binfac/N0)*exp(xset*tsol(j));
        end
        
        % Perturb the q posterior for the new event
        lampert = diag(lamtdiag(end, :));
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        qev(i+1, :) = qev(i+1, :)*lampert./sum(qev(i+1, :)*lampert);
        
        % Estimate the parameters of the rate and population
        xhatset{i} = zeros(1, length(tsol));
        xhatset2{i} = zeros(1, length(tsol));
        lamset{i} = zeros(1, length(tsol));
        thetaset{i} = zeros(1, length(tsol));
        Nset{i} = zeros(1, length(tsol));
        for j = 1:length(tsol)
            % Parameters, rate, 1/N(t) and N(t)
            lamset{i}(j) = qsol(j, :)*lamtdiag(j, :)';
            thetaset{i}(j) = (1/binfac)*qsol(j, :)*lamtdiag(j, :)';
            xhatset{i}(j) = qsol(j, :)*xset';
            Nset{i}(j) = binfac/(qsol(j, :)*lamtdiag(j, :)');
            
            % Other estimate from inverting lam
            xhatset2{i}(j) = (1/tsol(j))*log(lamset{i}(j)*N0/binfac);
        end
    end
    
    
    % Get full length of ODE solution data and assign appending vectors
    lenFull = sum(elemLen);
    stop = 0;
    qn = -ones(lenFull, m);
    tn = -ones(lenFull, 1);
    xhat = -ones(lenFull, 1);
    xhat2 = -ones(lenFull, 1);
    thetahat = -ones(lenFull, 1);
    lamhat = -ones(lenFull, 1);
    Nhat = -ones(lenFull, 1);
    
    % Append the cell based posterior and time data into a single structure
    for i = 1:nData
        % Loop calculates the start and end points along the array successively
        % and then assigns the appropriate cell element
        start = stop + 1;
        stop = elemLen(i) + start - 1;
        tn(start:stop) = tset{i};
        qn(start:stop, :) = qset{i};
        xhat(start:stop, :) = xhatset{i};
        xhat2(start:stop, :) = xhatset2{i};
        thetahat(start:stop, :) = thetaset{i};
        lamhat(start:stop, :) = lamset{i};
        Nhat(start:stop, :) = Nset{i};
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % Calculate true non-lineage dependent rate and population
    N = N0*exp(-x*tn);
    theta = 1./N;
    tnmin = min(tn);
    tnmax = max(tn);
    
%     L = 10;
%     F = 100;
%     [tx, Nhatx, lamhatx] = makeUniformSamples(tn, L, F, Nhat, lamhat);
%     Nx = N0*exp(-x*tx);
    
    % Calculate poisson rate for coalescent events %% <<<<<< check for repeats
    % in values due to duplicate tn values
    lam = zeros(size(tn));
    for i = 1:nData
        id = find(tn < tcoal(i+1) & tn >= tcoal(i));
        tni = tn(id);
        lam(id) = fac(i)*theta(id);
    end
    
    % Get the maximum and mean rate - by entering 0 one gets lam stats
    [lamstat, lam_m, tm] = coalMSE(tn, lam, zeros(size(lam)));
    lam_mean = lamstat.val(1);
    lam_max = max(lam_m);
    r = lam_max/lam_mean;
    C = lam_mean*log(r);
    disp(['Log of max to mean rate is ' num2str(log(r))]);
    
    % Plot the parameter,theta = 1/N, N estimates and rate
    if plotindiv
        figure;
        plot(tn, xhat, tn, x*ones(size(tn)), 'r');
        xlabel('time');
        ylabel('x in N(t) = N_0exp(-xt)');
        legend('xsny', 'x', 'location', 'best');
        title(['Parameter estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
        ylim([0.95*min([x; xhat]) 1.05*max([x; xhat])]);
        figure;
        plot(tn, thetahat, tn, theta, 'r');
        xlabel('time');
        ylabel('\theta');
        legend('thetasny', 'theta', 'location', 'best');
        title(['1/N(t) estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
        figure;
        semilogy(tn, lamhat, tn, lam, 'r');
        xlabel('time');
        ylabel('\lambda');
        legend('lamsny', 'lam', 'location', 'best');
        title(['Rate estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
        figure;
        loglog(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel('N(t)');
        legend('Nsny', 'N', 'location', 'best');
        title(['N(t) estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
    end
    
    % MSE ratio calculation and time scales and actual MSE
    [Nstat, em, tm, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
    R = 100*Nstat.val(3);
    [MSEstat, em, tm, ~] = coalMSE(tn, N, Nhat);
    
    % Remove zeros from null sets in running relative MSE and get length
    Rsqset = setData{1};
    setLen = setData{2};
    Rsqset = 100*Rsqset(Rsqset > 0);
    remSets = setLen(Rsqset <= 0);
    Rsqlen = 1:length(Rsqset);
    
    % Plot squared errors and MSE (perc ratios)
    if plotindiv
        figure;
        plot(Rsqlen, Rsqset, Rsqlen, R*ones(size(Rsqlen)), 'k');
        xlim([1 Rsqlen(end)]);
        xlabel('time');
        ylabel('% (Nhat/N - 1)^2');
        legend('cumulative on subsets', 'complete', 'location', 'best');
        title('View of the MSE ratio convergence');
    end
    
    % Store batch results
    MSEpercratio(ii) = R;
    MSE(ii) = MSEstat.val(3);
    param(ii) = x;
    paramhat(ii) = xhat(end);
    Cset(ii) = C;
    rset(ii) = r;
    disp('***************************************************************');
    disp(['The MSE perc ratio is ' num2str(R)]);
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
end

% Save important bits of data
if M > 1
    save(['batchExpdata2_' num2str(M)], 'MSEpercratio', 'param', 'paramhat',...
        'xmax', 'xmin', 'm', 'n', 'N0', 'fac', 'M');
    
    % Remove NaN values
    idrem1 = find(~isnan(MSEpercratio));
    idrem2 = find(~isnan(paramhat));
    idrem = find(~isnan(MSEpercratio) | ~isnan(paramhat));
    MSEpercratio = MSEpercratio(idrem);
    param = param(idrem);
    paramhat = paramhat(idrem);
    Cset = Cset(idrem);
    rset = rset(idrem);
    
    % Plot and histogram the percent MSE ratios and their stats
    Mperc = mean(MSEpercratio);
    Sperc = std(MSEpercratio);
    dom = 1:length(MSEpercratio);
    
    figure;
    plot(dom, MSEpercratio, dom, Mperc*ones(size(dom)), 'k', ...
        dom, (Mperc + 2*Sperc)*ones(size(dom)), 'r',...
        dom, (Mperc - 2*Sperc)*ones(size(dom)), 'r');
    xlabel('runs');
    ylabel('% (Nhat/N - 1)^2');
    title(['MSE performance for [n m mean %] = [' num2str(n) ' ' num2str(m) ' ' num2str(Mperc) ']']);
    
    figure;
    histfit(MSEpercratio, [], 'exp');
    xlabel('% (Nhat/N - 1)^2');
    ylabel('frequency');
    title('Fit of MSE performance to exponential');
    
    figure;
    eparam = param - paramhat;
    histfit(eparam, [], 'normal');
    xlabel('param - paramhat');
    ylabel('frequency');
    title('Fit of parameter error to normal');
end