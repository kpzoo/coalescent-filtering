% Perform Snyder filtering for for variable N(t) = xsin(wt) + xmax or the
% inverse depending on lamtype - problem with the inverse is NHPP sim

% Assumptions and Modifications
% - added repeatability for M runs
% - sinusoidal formulation for lam
% - the effective population varies with r parameters with m_i values
% - the data is k-coalescent times and the process self-exciting
% - modified to account for N being sinusoidal instead of theta

clc
close all
clear all

% Booleans to control functionality
M = 500;
plotindiv = 0;


% Parameters to store results
MSEpercratio = zeros(1, M);
param = zeros(1, M);
paramhat = zeros(1, M);

for ii = 1:M
    
    %% Generate the appropriate coalescent data for a defined N(t)
    
    % Set initial parameters: q0 is prior, N is true value, m is dimension
    % of parameter space, n is no. samples and m is total vector length
    m = 100;
    q0 = ones(1, m)/m;
    n = 100;
    nData = n-1;
    
    % Set population parameters for sinusoid
    xmax = n/1000000;
    xmin = 0.1*xmax;
    xset = linspace(xmin, xmax, m);
    x = randsample(xset, 1);
    w = 10;
    
    % Get coalescent binomial factors for lineages
    nset = n:-1:2;
    fac = nset.*(nset-1)/2;
    
    
    % % Simulate a homogeneous Poisson process with 50n points at max rate
    % Nhpp = 50*n;
    % lam_max = x + xmax;
    % thpp = exprnd(1/lam_max, Nhpp, 1);
    % thpp = cumsum(thpp);
    %
    % % Thin the homogeneous process to obtain the correct intensity
    % randNos = rand(Nhpp, 1);
    % rateNonHomo = x*sin(w*thpp) + xmax;
    % rateRatio = rateNonHomo/lam_max;
    % tcoal = thpp(randNos < rateRatio);
    % tcoal = [0 tcoal(1:nData)'];
    % tmrca = tcoal(end);
    % twait = diff(tcoal);
    
    % Generate the coalescent data with the appropriate rate
    lamtype = 1;
    [twait, tcoal] = getCoalData1Param(fac, x, nData, w, xmax, lamtype);
    tmrca = tcoal(end);
    twait2 = twait.*fac;
    tcoal2 = cumsum([0 twait2]);
    
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
        [tsol, qsol] = ode113(@(ts, y) odeSnyLinNonHomo1Param(ts, y, [], binfac, xset, w, lamtype),...
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
            % Calculate inhomogeneous function of choice
            switch(lamtype)
                case 1
                    % Sinusoidal theta(t)
                    lamtdiag(j, 1:m) = binfac*(xset*sin(w*tsol(j)) + xmax);
                case 2
                    % Sinusoidal N(t)
                    lamtdiag(j, 1:m) = binfac./(xset*sin(w*tsol(j)) + xmax);
            end
        end
        
        % Perturb the q posterior for the new event
        lampert = diag(lamtdiag(end, :));
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        qev(i+1, :) = qev(i+1, :)*lampert./sum(qev(i+1, :)*lampert);
        
        % Estimate the parameters of the rate and population
        xhatset{i} = zeros(1, length(tsol));
        lamset{i} = zeros(1, length(tsol));
        thetaset{i} = zeros(1, length(tsol));
        Nset{i} = zeros(1, length(tsol));
        for j = 1:length(tsol)
            % Parameters, rate and 1/N(t)
            lamset{i}(j) = qsol(j, :)*lamtdiag(j, :)';
            thetaset{i}(j) = (1/binfac)*qsol(j, :)*lamtdiag(j, :)';
            xhatset{i}(j) = qsol(j, :)*xset';
            Nset{i}(j) = binfac/(qsol(j, :)*lamtdiag(j, :)');
        end
    end
    
    
    % Get full length of ODE solution data and assign appending vectors
    lenFull = sum(elemLen);
    stop = 0;
    qn = -ones(lenFull, m);
    tn = -ones(lenFull, 1);
    xhat = -ones(lenFull, 1);
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
        thetahat(start:stop, :) = thetaset{i};
        lamhat(start:stop, :) = lamset{i};
        Nhat(start:stop, :) = Nset{i};
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % Calculate true non-lineage dependent rate and population
    switch(lamtype)
        case 1
            % Sinusoidal theta(t)
            theta = x*sin(w*tn) + xmax;
            N = 1./(x*sin(w*tn) + xmax);
        case 2
            % Sinusoidal N(t)
            N = x*sin(w*tn) + xmax;
            theta = 1./(x*sin(w*tn) + xmax);
    end
    tnmin = min(tn);
    tnmax = max(tn);
    
    lam = zeros(size(tn));
    % Calculate poisson rate for coalescent events %% <<<<<< check for repeats
    % in values due to duplicate tn values
    for i = 1:nData
        id = find(tn < tcoal(i+1) & tn >= tcoal(i));
        tni = tn(id);
        switch(lamtype)
            case 1
                % Sinusoidal theta(t)
                lam(id) = fac(i)*(x*sin(w*tni) + xmax);
            case 2
                % Sinusoidal N(t)
                lam(id) = fac(i)/(x*sin(w*tni) + xmax);
        end
    end
    
    % Get the maximum and mean rate - by entering 0 one gets lam stats
    [lamstat, lam_m, tm] = coalMSE(tn, lam, zeros(size(lam)));
    lam_mean = lamstat.val(1);
    lam_max = max(lam_m);
    r = lam_max/lam_mean;
    disp(['Log of max to mean rate is ' num2str(log(r))]);
    
    % Plot the parameter and theta = 1/N estimates and rate
    if plotindiv
        figure;
        plot(tn, xhat, tn, x*ones(size(tn)), 'r');
        xlabel('time');
        ylabel('x in xsin(wt)+xmax');
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
        semilogy(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel('N(t)');
        legend('Nsny', 'N', 'location', 'best');
        title(['N(t) estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
    end
    
    % MSE ratio calculation and time scales
    [Nstat, em, tm, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
    R = 100*Nstat.val(3);
    
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
    param(ii) = x;
    paramhat(ii) = xhat(end);
    disp('***************************************************************');
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
end

% Save important bits of data
save(['batchdata_' num2str(M)], 'MSEpercratio', 'param', 'paramhat',...
    'xmax', 'xmin', 'm', 'n', 'w', 'fac', 'M');

% Plot and histogram the percent MSE ratios and their stats
Mperc = mean(MSEpercratio);
Sperc = std(MSEpercratio);
dom = 1:M;

figure;
plot(dom, MSEpercratio, dom, Mperc*ones(size(dom)), 'k', ...
    dom, (Mperc + 2*Sperc)*ones(size(dom)), 'r',...
    dom, (Mperc - 2*Sperc)*ones(size(dom)), 'r');
xlabel('runs');
ylabel('% (Nhat/N - 1)^2');
title(['MSE performance for [n m mean %] = [' num2str(n) ' ' num2str(m) ' ' num2str(Mperc) ']']);

figure;
histfit(MSEpercratio, 25, 'exp');
xlabel('% (Nhat/N - 1)^2');
ylabel('frequency');
title('Fit of MSE performance to exponential');