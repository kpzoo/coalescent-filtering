% Perform Snyder filtering for constant N only using standard measures for
% comparison and placement in technical report

% Assumptions and Modifications
% - removed the scaling that gave Nev2 etc - no time scaling
% - the effective population is constant and can take m random values
% - the data is k-coalescent times and the process self-exciting
% - added option to run M times to get comparisons

clc
close all
clear all

% Booleans to control functionality
M = 1;
plotindiv = 1;

% Parameters to store results
Nval = zeros(1, M);
Npar = zeros(1, M);
Nnopar = zeros(1, M);

for ii = 1:M
    
    %% Filter the simulated event times from the appropriate Poisson process
    
    % Set initial parameters: q0 is prior, N is true value, m is dimension of N
    % space and n is no. samples
    m = 100;
    q0 = ones(1, m)/m;
    n = 100;
    nData = n-1;
    
    % Get the instance of random variable N and initial rate set vector without
    % the self excitation, lvec
    %Nset = 50*(1:m) + 100;
    Nset = n*linspace(100, 1000, m);
    N = randsample(Nset, 1);
    lvec = 1./Nset;
    
    % Get rates of coalescent and the waiting times, twait as well as
    % cumulative times and the MRCA time
    nset = n:-1:2;
    fac = nset.*(nset-1)/2;
    rates = fac/N;
    twait = exprnd(1./rates);
    tcoal = cumsum([0 twait]);
    tmrca = tcoal(end);
    
    % Posterior vectors on events, qnoev is for non-updated event q's
    qev = zeros(nData+1, m);
    qnoev = qev;
    qev(1, :) = q0;
    tev = zeros(nData, 1);
    
    % Cell to save output of ODE solver and set options
    qset = cell(1, 1);
    tset = cell(1, 1);
    Nhatset = cell(1, 1);
    options = odeset('NonNegative', 1:m);
    elemLen = zeros(1, nData);
    Nev = zeros(1, nData);
    
    % Run Snyder filter across the time series
    for i = 1:nData
        % Obtain rate matrix, lam with account for scaling
        lam = diag(fac(i)*lvec);
        
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam),...
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
        
        % Perturb the q posterior for the new event
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        qev(i+1, :) = qev(i+1, :)*lam./sum(qev(i+1, :)*lam);
        
        % Estimate the population parameter N from the posterior
        Nhatset{i} = fac(i)./(qsol*diag(lam));
        Nev(i) = fac(i)./(qev(i+1, :)*diag(lam));
    end
    
    
    % Get full length of ODE solution data and assign appending vectors
    lenFull = sum(elemLen);
    stop = 0;
    qn = -ones(lenFull, m);
    tn = -ones(lenFull, 1);
    Nhat = -ones(lenFull, 1);
    
    % Append the cell based posterior and time data into a single structure
    for i = 1:nData
        % Loop calculates the start and end points along the array successively
        % and then assigns the appropriate cell element
        start = stop + 1;
        stop = elemLen(i) + start - 1;
        tn(start:stop) = tset{i};
        qn(start:stop, :) = qset{i};
        Nhat(start:stop, :) = Nhatset{i};
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % MSE calculation and time scales
    [statSet, em, tm] = coalMSE(tn, N, Nhat);
    tnmin = min(tn);
    tnmax = max(tn);
    
    if plotindiv
        % Plot the estimate ratio of N
        figure;
        plot(tn, Nhat/N, [tnmin; tnmax], [1; 1], 'k');
        ylim([0.95 1.05]);
        xlabel('time');
        ylabel('Nhat/N estimate');
        title(['Snyder estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
    end
    
    % Check against classical MMSE estimators = suff stat/(ndata - 2)
    T1 = sum(fac.*twait);
    Nclass = T1/(nData - 2);
    Nsny = Nhat(end);
    
    % Display comparisons
    disp(['N is ' num2str(N)]);
    disp(['Nhat(end) is ' num2str(Nsny)]);
    disp(['Nclass is ' num2str(Nclass)]);
    
    % Get squared error ratios
    RsqErr = [(1 - Nsny/N)^2 (1 - Nclass/N)^2];
    disp(['Squared error ratios are: ' num2str(RsqErr)]);
    
    % Get cumulative classical estimators
    t1 = cumsum(fac.*twait);
    rr = (1:nData) - 2;
    Nc = t1./rr;
    
    % Remove the division by zero and negative values and make Nevs consistent
    Nc = Nc(3:end);
    Nevf = Nev(3:end);
    newlen = length(Nc);
    
    % Get MSE ratio of cumulative classical and Snyder estimators
    mseCum = zeros(2, newlen);
    mseCum(1, :) = (1 - Nevf/N).^2;
    mseCum(2, :) = (1 - Nc/N).^2;
    
    if plotindiv
        % Plot estimators across events or data samples
        figure;
        semilogy(1:newlen, mseCum(1, :), 'b', 1:newlen, mseCum(2, :), 'r');
        xlabel('events');
        ylabel('relative MSE');
        legend('param', 'non-param', 'location', 'best');
        title(['Relative estimator performance for [n m] = [' num2str(n) ' ' num2str(m) ']']);
    end
    
    % Save key results
    Nval(ii) = N;
    Npar(ii) = Nsny;
    Nnopar(ii) = Nclass;
    disp(['Finished iteration: ' num2str(ii) ' of ' num2str(M)]);
    
end

%% Final analysis of runs
if M > 1
    save(['comp' num2str(n) num2str(m)]);
    
    % Square errors (percentage and relative)
    Spar = 100*((Nval - Npar).^2)./(Nval.^2);
    Snopar = 100*((Nval - Nnopar).^2)./(Nval.^2);
    Meanpar = mean(Spar);
    Meannopar = mean(Snopar);
    Varpar = var(Spar);
    Varnopar = var(Snopar);
    L = 1:M;
    
    % Plot square errors
    figure;
    stairs(L', [Snopar' Spar']);
    xlabel('runs');
    ylabel('relative percent MSE');
    legend('non-param', 'snyder', 'location', 'best');
    title(['Non-parametric and Snyder performance for n =' num2str(n)]);
    
    % Error difference between estimators
    figure;
    err = Snopar - Spar;
    stairs(err);
    ylabel('relative percent MSE');
    xlabel('runs');
    title(['Difference of square errors for [n m] = [' num2str(n) ' ' num2str(m) ']']);
end