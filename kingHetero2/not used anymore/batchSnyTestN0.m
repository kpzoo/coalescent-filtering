% Perform Snyder filtering for constant N(t) = N(0) = x [only for this N]

% Assumptions and Modifications
% - code only works for constant population <----------------------
% - added repeatability for M runs
% - tested with reversal of data and no scaling
% - the data is k-coalescent times and the process self-exciting
% - includes timing of simulations
% - removed thetahat so it is not calculated
% - removed capacity calculations

clc
close all
clear all

%profile on

% Booleans to control functionality
tic;
M = 1000;
plotindiv = 0;
reversal = 1;
tIter = zeros(1, M);

% Parameters to store results
MSEpercratio = zeros(1, M);
param = zeros(1, M);
paramhat = zeros(1, M);
if reversal
    MSEpercratio2 = zeros(1, M);
    param2 = zeros(1, M);
    paramhat2 = zeros(1, M);
end

for ii = 1:M
    
    %% Generate the appropriate coalescent data for a defined N(t)
    
    % Set initial parameters: q0 is prior, N is true value, m is dimension
    % of parameter space, n is no. samples and m is total vector length
    m = 100;
    q0 = ones(1, m)/m;
    n = 100;
    nData = n-1;
    
    % Set population parameters for exponential with N0 >> n
    xmax = 1000*n;
    xmin = 100*n;
    xset = linspace(xmin, xmax, m);
    x = randsample(xset, 1);
    if reversal
        xset2 = xset;
        x2 = x;
    end
    
    % Get coalescent binomial factors for lineages
    nset = n:-1:2;
    fac = nset.*(nset-1)/2;
    if reversal
        fac2 = fac(end:-1:1);
    end
    
    % Get set of lam (rate matrix) diagonals
    lamDiag = zeros(nData, m);
    for i = 1:nData
        lamDiag(i, :) = fac(i)./xset;
    end
    
    % Simulate coalescent times from appropriate exponential distribution
    rates = fac/x;
    twait = exprnd(1./rates);
    tcoal = cumsum([0 twait]);
    
    % Reorder the wait times and lamDiag for reversal of same data
    if reversal
        twait2 = twait(end:-1:1);
        tcoal2 = cumsum([0 twait2]);
        rates2 = rates(end:-1:1);
        lamDiag2 = lamDiag(end:-1:1, :);
    end
    
    %% Filter the simulated event times from the appropriate Poisson process
    
    % Posterior vectors on events, qnoev is for non-updated event q's
    qev = zeros(nData+1, m);
    qnoev = qev;
    qev(1, :) = q0;
    tev = zeros(nData, 1);
    if reversal
        qev2 = zeros(nData+1, m);
        qnoev2 = qev2;
        qev2(1, :) = q0;
        tev2 = zeros(nData, 1);
    end
    
    % Options setting for ODE solver
    options = odeset('NonNegative', 1:m);
    
    % Variables to save output of ODE solver
    qset = cell(1, 1);
    tset = cell(1, 1);
    xhatset = cell(1, 1);
    lamhatset = cell(1, 1);
    elemLen = zeros(1, nData);
    xev = zeros(1, nData);
    
    % Similar ODE outputs for test, which can be reversal or other specs+
    if reversal
        qset2 = cell(1, 1);
        tset2 = cell(1, 1);
        xhatset2 = cell(1, 1);
        lamhatset2 = cell(1, 1);
        elemLen2 = zeros(1, nData);
        xev2 = zeros(1, nData);
    end
    
    % Run Snyder filter across the time series
    for i = 1:nData
        % Obtain rate matrices and binomial factors
        binfac = fac(i);
        lam = diag(lamDiag(i, :));
        if reversal
            binfac2 = fac2(i);
            lam2 = diag(lamDiag2(i, :));
        end
        
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam),...
            [tcoal(i) tcoal(i+1)], qev(i, :)', options);
        if reversal
            [tsol2, qsol2] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam2),...
                [tcoal2(i) tcoal2(i+1)], qev2(i, :)', options);
        end
        
        % Normalise the posterior probabilities
        for j = 1:size(qsol, 1)
            qsol(j, :) = qsol(j, :)/(sum(qsol(j, :)));
        end
        if reversal
            for j = 1:size(qsol2, 1)
                qsol2(j, :) = qsol2(j, :)/(sum(qsol2(j, :)));
            end
        end
        
        % Assign the output values of time and posteriors
        qnoev(i, :) = qev(i, :);
        qset{i} = qsol;
        tset{i} = tsol;
        elemLen(i) = length(tsol);
        if reversal
            qnoev(i, :) = qev(i, :);
            qset2{i} = qsol2;
            tset2{i} = tsol2;
            elemLen2(i) = length(tsol2);
        end
        
        % Perturb the q posterior for the new event
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        qev(i+1, :) = qev(i+1, :)*lam./sum(qev(i+1, :)*lam);
        if reversal
            qev2(i+1, :) = qsol2(end, :);
            tev2(i+1) = tsol2(end);
            qev2(i+1, :) = qev2(i+1, :)*lam2./sum(qev2(i+1, :)*lam2);
        end
        
        % Estimate the parameters of the rate and population
        xhatset{i} = zeros(1, length(tsol));
        lamhatset{i} = zeros(1, length(tsol));
        for j = 1:length(tsol)
            % Parameters, rate and N(t)
            lamhatset{i}(j) = qsol(j, :)*diag(lam);
            xhatset{i}(j) = qsol(j, :)*xset';
        end
        
        if reversal
            xhatset2{i} = zeros(1, length(tsol2));
            lamhatset2{i} = zeros(1, length(tsol2));
            for j = 1:length(tsol2)
                % Parameters, rate, 1/N(t) and N(t)
                lamhatset2{i}(j) = qsol2(j, :)*diag(lam2);
                xhatset2{i}(j) = qsol2(j, :)*xset2';
            end     
        end
        
    end
    
    
    % Get full length of ODE solution data and assign appending vectors
    lenFull = sum(elemLen);
    stop = 0;
    qn = -ones(lenFull, m);
    tn = -ones(lenFull, 1);
    xhat = -ones(lenFull, 1);
    lamhat = -ones(lenFull, 1);
    if reversal
        lenFull2 = sum(elemLen2);
        stop2 = 0;
        qn2 = -ones(lenFull2, m);
        tn2 = -ones(lenFull2, 1);
        xhat2 = -ones(lenFull2, 1);
        lamhat2 = -ones(lenFull2, 1);
    end
    
    % Append the cell based posterior and time data into a single structure
    for i = 1:nData
        % Loop calculates the start and end points along the array successively
        % and then assigns the appropriate cell element
        start = stop + 1;
        stop = elemLen(i) + start - 1;
        tn(start:stop) = tset{i};
        qn(start:stop, :) = qset{i};
        xhat(start:stop, :) = xhatset{i};
        lamhat(start:stop, :) = lamhatset{i};
        if reversal
            start2 = stop2 + 1;
            stop2 = elemLen2(i) + start2 - 1;
            tn2(start2:stop2) = tset2{i};
            qn2(start2:stop2, :) = qset2{i};
            xhat2(start2:stop2, :) = xhatset2{i};
            lamhat2(start2:stop2, :) = lamhatset2{i};
        end
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % Calculate true population
    N = x*ones(size(tn));
    Nhat = xhat;
    tnmin = min(tn);
    tnmax = max(tn);
    if reversal
        N2 = x2*ones(size(tn2));
        Nhat2 = xhat2;
        tnmin2 = min(tn2);
        tnmax2 = max(tn2);
        % Overall time range
        trange = [min([tnmin, tnmin2]), max([tnmax, tnmax2])];
    end
    
    % Calculate poisson rate for coalescent events %% <<<<<< check for repeats
    % in values due to duplicate tn values
    lamn = zeros(size(tn));
    for i = 1:nData
        id = find(tn < tcoal(i+1) & tn >= tcoal(i));
        tni = tn(id);
        lamn(id) = fac(i)/x;
    end
    if reversal
        lamn2 = zeros(size(tn2));
        for i = 1:nData
            id2 = find(tn2 < tcoal2(i+1) & tn2 >= tcoal2(i));
            tni2 = tn2(id2);
            lamn2(id2) = fac2(i)/x2;
        end
    end
    
    % Twenty iterations of the posterior from prior
    nPoster = 20;
    finidpost = size(qn, 1);
    idpost = round(linspace(1, finidpost, nPoster));
    qnp = qn(idpost, :);
    if reversal
        finidpost2 = size(qn2, 1);
        idpost2 = round(linspace(1, finidpost2, nPoster));
        qnp2 = qn2(idpost2, :);
    end
    
    % MSE ratio calculation and time scales and actual MSE
    [Nstat, ~, ~, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
    R = 100*Nstat.val(3);
    
    if reversal
        [Nstat2, ~, ~, setData2] = coalMSE(tn2, Nhat2./N2, ones(size(tn2)));
        R2 = 100*Nstat2.val(3);
    end
    
    % Remove zeros from null sets in running relative MSE and get length
    Rsqset = setData{1};
    setLen = setData{2};
    Rsqset = 100*Rsqset(Rsqset > 0);
    remSets = setLen(Rsqset <= 0);
    Rsqlen = 1:length(Rsqset);
    if reversal
        Rsqset2 = setData2{1};
        setLen2 = setData2{2};
        Rsqset2 = 100*Rsqset2(Rsqset2 > 0);
        remSets2 = setLen2(Rsqset2 <= 0);
        Rsqlen2 = 1:length(Rsqset2);
    end
    
    %% Perform plotting and visualisation
    
    % Plot the individual non-reversed results
    if plotindiv && ~reversal
        % The complete coalescent rate
        figure;
        semilogy(tn, lamhat, tn, lamn, 'r');
        xlabel('time');
        ylabel('\lambda');
        legend('lamsny', 'lam', 'location', 'best');
        title(['Rate estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
        
        % The population function
        figure;
        plot(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel('N(t)');
        legend('Nsny', 'N_0', 'location', 'best');
        title(['N(t) estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim([tnmin tnmax]);
        
        if reversal
            figure;
            semilogy(tn2, lamhat2, tn2, lamn2, 'r');
            xlabel('time');
            ylabel('\lambda (reversed)');
            legend('lamsny', 'lam', 'location', 'best');
            title(['Reverse rate estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
            xlim([tnmin2 tnmax2]);
        end
        
        % Plot the posterior evolution
        figure;
        plot(xset, qnp(1, :), 'b', xset, qnp(end, :), 'k',...
            xset, qnp(2:end-1, :), 'g--');
        hold on
        plot([x x], [0 max(max(qnp))], 'k');
        hold off
        xlabel('N_0');
        ylabel('P(N_0|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors with N_0 = ' num2str(x)]);
        
        % Plot squared errors and MSE (perc ratios)
        figure;
        plot(Rsqlen, Rsqset, Rsqlen, R*ones(size(Rsqlen)), 'k');
        xlim([1 Rsqlen(end)]);
        xlabel('subsets');
        ylabel('% (Nhat/N_0 - 1)^2');
        legend('cumulative on subsets', 'complete', 'location', 'best');
        title('View of the MSE ratio convergence');
       
    end
    
    % Specific reversal plots - look at final distributions
    if reversal && plotindiv
        figure;
        plot(xset, qn(end, :), 'bo', xset2, qn2(end, :), 'ro');
        hold on
        plot([Nhat(end) Nhat(end)], [0 max(qn(end, :))], 'g');
        plot([Nhat2(end) Nhat2(end)], [0 max(qn2(end, :))], 'm');
        hold off
        xlabel('N_0');
        ylabel('P(N_0|data)');
        legend('normal', 'reversal', 'Nhat norm', 'Nhat rev', 'location', 'best');
        title('Posterior distributions for reversed and normal data');
        
        % Compare posterior distributions and final estimate
        figure;
        subplot(2, 1, 1);
        plot(xset, qnp(1, :), 'b', xset, qnp(end, :), 'k',...
            xset, qnp(2:end-1, :), 'g--');
        hold on
        plot([x x], [0 max(max(qnp))], 'k');
        hold off
        xlabel('N_0');
        ylabel('P(N_0|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: [N_0 Nhat] = ' [num2str(x) ', ' num2str(xhat(end))]]);
        subplot(2, 1, 2);
        plot(xset2, qnp2(1, :), 'b', xset2, qnp2(end, :), 'k',...
            xset2, qnp2(2:end-1, :), 'g--');
        hold on
        plot([x2 x2], [0 max(max(qnp2))], 'k');
        hold off
        xlabel('N_0');
        ylabel('P(N_0|data) reversed');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of reversed posteriors: [N_0 Nhat] = ' [num2str(x2) ', ' num2str(xhat2(end))]]);
       
        % Compare square error convergence
        figure;
        plot(Rsqlen, Rsqset, Rsqlen2, Rsqset2);
        xlabel('subsets');
        ylabel('% (Nhat/N_0 - 1)^2');
        legend('normal', 'reversed', 'location', 'best');
        title('Comparison of the MSE ratio convergence');
        
        % Population functions
        figure;
        plot(tn, Nhat, tn2, Nhat2, trange, [x x], 'k');
        xlabel('time');
        ylabel('N(t)');
        legend('normal', 'reversed', 'true', 'location', 'best');
        title(['N_0 estimates for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim(trange);
        
        % The complete coalescent rate
        figure;
        semilogy(tn, lamhat, tn, lamn, tn2, lamhat2, tn2, lamn2);
        xlabel('time');
        ylabel('\lambda');
        legend('est normal', 'true normal', 'est reversed', 'true reversed', 'location', 'best');
        title(['Coalescent rate estimates for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        xlim(trange);
    end
    
    %% Store batch results
    
    % Main variables of interest
    MSEpercratio(ii) = R;
    param(ii) = x;
    paramhat(ii) = xhat(end);
    if reversal
        MSEpercratio2(ii) = R2;
        param2(ii) = x2;
        paramhat2(ii) = xhat2(end);
    end
    
    % Display progress and clock iteration time
    tIter(ii) = toc;
    disp('***************************************************************');
    disp(['The MSE perc ratio is ' num2str(R)]);
    if reversal
        disp(['The MSE perc reverse ratio is ' num2str(R2)]);
    end
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
end

% Clock time set to minutes
tIter = tIter/60;
tavgIter = mean(diff(tIter));
disp(['Total time for ' num2str(M) ' iterations is ' num2str(tIter(end))]);
disp(['Average time with [n m] = ' [num2str(n) ' ' num2str(m)] ' is ' num2str(tavgIter)]);

%profile viewer

% Save the results and plot batch data
if M > 1
    % Save data
    name = ['NconstRev' num2str(M)];
    save(name, 'MSEpercratio', 'param', 'paramhat', 'tIter',...
        'xmax', 'xmin', 'm', 'n', 'x', 'fac', 'M');
    if reversal
        % Reversal data to append
        save(name, 'MSEpercratio2', 'param2', 'paramhat2', 'fac2', '-append');
    end
    
    % Fitting of parameter estimates to normal
    figure;
    eparam = param - paramhat;
    histfit(eparam, [], 'normal');
    xlabel('param - paramhat');
    ylabel('frequency');
    title('Fit of parameter error to normal');
    if reversal
        figure;
        eparam2 = param2 - paramhat2;
        histfit(eparam2, [], 'normal');
        xlabel('param - paramhat');
        ylabel('frequency');
        title('Fit of parameter error to normal for reversed data');
    end
    
    % More sensible relative MSE comes from final value of parameters when
    % N_0 = x (constant case only) <---------------------
    relmse = 100*esq/(x^2);
    mr = mean(relmse)*ones(1, M);
    sr = std(relmse)*ones(1, M);
    
    figure;
    plot(1:M, relmse, 1:M, mr, 'k', 1:M, mr - 2*sr, 'r', 1:M, mr + 2*sr, 'r')
    xlabel('runs')
    ylabel('% (1 - Nhat/N_0)^2')
    title(['Relative last value MSE across runs, mean = ' num2str(mr(1))]);
    legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    
    if reversal
        % Calculate MSE in parameters
        esq = eparam.^2;
        esq2 = eparam2.^2;
        MSE1 = mean(esq);
        MSE2 = mean(esq);
        
        % Plot the correlation comparison
        figure;
        plot(eparam, eparam2, 'ro', eparam, eparam, 'k');
        xlabel('N_0 - E[N_0 | data]');
        ylabel('N_0 - E[N_0 | data], reverse time');
        legend('error comparison', 'y = x line', 'location', 'best');
        title('Perfect correlation between normal and reverse time estimates');
    end
end