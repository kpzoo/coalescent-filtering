% Perform Snyder filtering for variable N(t) = x1exp(-x2t) (reverse time)
% where x1 = N0 and x2 = r in standard notation

% Assumptions and Modifications
% - removed theta and capacity and MSEperc to MSE of parameters across runs
% - added repeatability for M runs
% - exponential formulation for lam but with 2 params
% - used time rescaling from Palacios 2012 to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clear all

% Booleans to control functionality
tic;
M = 1;
plotindiv = 1;
tIter = zeros(1, M);

%% Initialise and define key parameters and storage variables

% Set number of parameters, dimensions (mi) abd data length
mi = [100 100];
numRV = length(mi);
n = 200;
nData = n-1;

% Set space, prior q0 and total dimension m
minspace = [100*n 0.1];
maxspace = [1000*n 10];
m = prod(mi);
q0 = ones(1, m)/m;

% Get the instance of the random variables for rate function, lam
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minspace(i), maxspace(i), mi(i));
    x(i) = randsample(xset{i}, 1);
end

% Get coalescent binomial factors for lineages and change notation
nset = n:-1:2;
fac = nset.*(nset-1)/2;
N0 = x(1);
r = x(2);

% Get (x1, x2) combined into a single array via kronecker products
rset = xset{2};
N0set = xset{1};
space = kron(N0set, rset);
    
% Parameters to store results
MSEpercratio = zeros(1, M);
param = zeros(numRV, M);
paramhat = zeros(numRV, M);

for ii = 1:M
    
    %% Generate the appropriate coalescent data for a defined N(t)

    % Simulate coalescent times for exponential using time rescaling
    tcoal = zeros(1, n);
    for i = 1:n-1
        % Get coalescent time for i+1 -> i lineages
        tcoal(i+1) = log(exp(r*tcoal(i)) + r*N0*exprnd(1)/fac(i))/r;
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
    lamhatset = cell(1, 1);
    Nhatset = cell(1, 1);
    options = odeset('NonNegative', 1:m);
    elemLen = zeros(1, nData);
    xev = zeros(1, nData);
    
    % Diagonal of non-time dependent but event dependent lambda
    lamnoTdiag = cell(1, nData);
    
    % Run Snyder filter across the time series
    for i = 1:nData
        % Obtain approproate binomial factor
        binfac = fac(i);
        
        % Calculate non-time dependent part of lambda
        lamnoTdiag{i} = repmat(binfac./xset{1}, 1, mi(2));
        
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnyExp2Param(ts, y, [], lamnoTdiag{i}, xset, mi),...
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
            for kk = 1:mi(2)
                % Index the start and end points based on mi(2)
                id1 = 1 + (kk - 1)*mi(1);
                id2 = kk*mi(1);
                lamtdiag(j, id1:id2) = lamnoTdiag{i}(id1:id2)*exp(xset{2}(kk)*tsol(j));
            end
        end
        
        % Perturb the q posterior for the new event
        lampert = diag(lamtdiag(end, :));
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        qev(i+1, :) = qev(i+1, :)*lampert./sum(qev(i+1, :)*lampert);
        
        % Estimate the rate and population directly
        lamhatset{i} = zeros(1, length(tsol));
        Nhatset{i} = zeros(1, length(tsol));
        for j = 1:length(tsol)
            % Parameters, rate and N(t)
            lamhatset{i}(j) = qsol(j, :)*lamtdiag(j, :)';
            Nhatset{i}(j) = binfac/(qsol(j, :)*lamtdiag(j, :)');
        end
    end
    
    
    % Get full length of ODE solution data and assign appending vectors
    lenFull = sum(elemLen);
    stop = 0;
    qn = -ones(lenFull, m);
    tn = -ones(lenFull, 1);
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
        lamhat(start:stop, :) = lamhatset{i};
        Nhat(start:stop, :) = Nhatset{i};
    end
    
    % Define cell to hold probabilities, phat for different RVs
    phat = cell(1);
    tnlen = length(tn);
    for kk = 1:numRV
        phat{kk} = zeros(tnlen, length(xset{kk}));
    end
    
    % Marginalise the posterior to obtain parameter estimates
    for j = 1:tnlen
        % Reshape the posterior in a form so that sums on each dimension
        % give marginal probabilities at each time
        qtemp = reshape(qn(j, :), mi);
        for kk = 1:numRV
            qsum = sum(qtemp, numRV - kk + 1);
            phat{kk}(j, :) = qsum;
        end
    end
    
    % Obtain parameter estimates
    xhat = zeros(tnlen, numRV);
    for kk = 1:numRV
        xhat(:, kk) = phat{kk}*xset{kk}';
    end
    
    % Obtain parameter estimates E[xi|data] and the std deviations -
    % separate vectors as orders of magnitude different
    N0hat = phat{1}*xset{1}';
    rhat = phat{2}*xset{2}';
    % Std deviation of conditional mean
    N0hatvar = phat{1}*(xset{1}'.^2) - N0hat.^2;
    rhatvar = phat{2}*(xset{2}'.^2) - rhat.^2;
    N0hatstd = sqrt(N0hatvar);
    rhatstd = sqrt(rhatvar);
    % Calculate upper and lower std deviation bounds mean +/- 2*std
    N0ubstd = N0hat + 2*N0hatstd;
    N0lbstd = N0hat - 2*N0hatstd;
    rubstd = rhat + 2*rhatstd;
    rlbstd = rhat - 2*rhatstd;
    
    %% Convert data in MSE estimates and plot comparisons
    
    % Calculate true non-lineage dependent rate and population
    N = N0*exp(-r*tn);
    tnmin = min(tn);
    tnmax = max(tn);
    
    % Calculate poisson rate for coalescent events %% <<<<<< check for repeats
    % in values due to duplicate tn values
    lam = zeros(size(tn));
    for i = 1:nData
        id = find(tn < tcoal(i+1) & tn >= tcoal(i));
        tni = tn(id);
        lam(id) = fac(i)./N(id);
    end
    
    % Twenty iterations of the posterior from prior
    nPoster = 20;
    qnp = cell(1, numRV);
    for kk = 1:numRV
        finidpost = tnlen;
        idpost = round(linspace(1, finidpost, nPoster));
        % Posteriors spaced as specified by idpost
        qnp{kk} = phat{kk}(idpost, :);
    end
    
    % Plot the N estimates and coalescent rate lambda
    if plotindiv
        % The coalescent rate which includes lineages
        figure;
        semilogy(tn, lamhat, tn, lam, 'r');
        xlabel('time');
        ylabel('\lambda(N_0, r, t)');
        legend('lamsny', 'lam', 'location', 'best');
        title(['Rate estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
        % The population function
        figure;
        semilogy(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel('N(t) = N_0e^{-rt}');
        legend('Nsny', 'N', 'location', 'best');
        title(['N(t) estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
        % The parameters estimated
        figure;
        subplot(2, 1, 1);
        plot(tn, x(1)*ones(size(tn)), 'r', tn, N0hat);
        xlabel('time');
        ylabel('N_0 in N(t) = N_0e^{-rt}');
        title(['Estimate of N_0 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        ylim([0.95*min([x(1); N0hat]) 1.05*max([x(1); N0hat])]);
        xlim([tnmin tnmax]);
        subplot(2, 1, 2);
        plot(tn, x(2)*ones(size(tn)), 'r', tn, rhat);
        xlabel('time');
        ylabel('r in N(t) = N_0e^{-rt}');
        title(['Estimate of r for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        ylim([0.95*min([x(2); rhat]) 1.05*max([x(2); rhat])]);
        xlim([tnmin tnmax]);
        
        % The parameters estimated and the std deviation bounds
        figure;
        subplot(2, 1, 1);
        plot(tn, x(1)*ones(size(tn)), 'r', tn, N0hat, 'k', tn, N0lbstd, 'b', tn, N0ubstd, 'b');
        xlabel('time');
        ylabel('N_0 in N(t) = N_0e^{-rt}');
        title(['Estimate of N_0 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        %ylim([0.95*min([x(1); N0hat]) 1.05*max([x(1); N0hat])]);
        legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
        xlim([tnmin tnmax]);
        subplot(2, 1, 2);
        plot(tn, x(2)*ones(size(tn)), 'r', tn, rhat, 'k', tn, rlbstd, 'b', tn, rubstd, 'b');
        xlabel('time');
        ylabel('r in N(t) = N_0e^{-rt}');
        title(['Estimate of r for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        %ylim([0.95*min([x(2); rhat]) 1.05*max([x(2); rhat])]);
        legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
        xlim([tnmin tnmax]);
        
        % Posteriors and final estimates
        figure;
        subplot(2, 1, 1);
        plot(xset{1}, qnp{1}(1, :), 'b', xset{1}, qnp{1}(end, :), 'k',...
            xset{1}, qnp{1}(2:end-1, :), 'g--');
        hold on
        plot([x(1) x(1)], [0 max(max(qnp{1}))], 'k');
        hold off
        xlabel('N_0');
        ylabel('P(N_0|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: [N_0 N_0hat] = ' [num2str(x(1)) ', ' num2str(N0hat(end))]]);
        subplot(2, 1, 2);
        plot(xset{2}, qnp{2}(1, :), 'b', xset{2}, qnp{2}(end, :), 'k',...
            xset{2}, qnp{2}(2:end-1, :), 'g--');
        hold on
        plot([x(2) x(2)], [0 max(max(qnp{2}))], 'k');
        hold off
        xlabel('r');
        ylabel('P(r|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: [r rhat] = ' [num2str(x(2)) ', ' num2str(rhat(end))]]);
    end
    
    % MSE ratio calculation as a percent
    [Nstat, ~, ~, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
    R = 100*Nstat.val(3);
    
    % Store batch results with param and paramhat having rows for numRVs
    MSEpercratio(ii) = R;
    param(:, ii) = x';
    paramhat(:, ii) = xhat(end, :)';
    disp('***************************************************************');
    disp(['The MSE perc ratio is ' num2str(R)]);
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
    tIter(ii) = toc;
end

% Clock time set to minutes
tIter = tIter/60;
disp(['Total time for ' num2str(M) ' iterations is ' num2str(tIter(end))]);
if M > 1
    tavgIter = mean(diff(tIter));
    disp(['Average time with [n m] = ' [num2str(n) ' ' num2str(m)] ' is ' num2str(tavgIter)]);
end

% Save important bits of data
if M > 1
    save(['2paramexp_' num2str(M)], 'MSEpercratio', 'param', 'paramhat',...
        'space', 'xset', 'm', 'n', 'mi', 'x', 'fac', 'M', 'tIter');
    
    % Parameter estimate errors in multivariate form
    eparam = param - paramhat;
    esq = eparam.^2;
    
    % Get specific estimates for N0 and r and their statistics
    relmse_N0 = 100*esq(1, :)/(x(1)^2);
    relmse_r = 100*esq(2, :)/(x(2)^2);
    mN0 = mean(relmse_N0)*ones(1, M);
    sN0 = std(relmse_N0)*ones(1, M);
    mr = mean(relmse_r)*ones(1, M);
    sr = std(relmse_r)*ones(1, M);
    
    % Plot the behaviour across runs
    figure;
    subplot(2, 1, 1);
    plot(1:M, relmse_N0, 1:M, mN0, 'k', 1:M, mN0 - 2*sN0, 'r', 1:M, mN0 + 2*sN0, 'r')
    xlabel('runs')
    ylabel('% (1 - Nhat/N_0)^2')
    title(['Relative MSE(end) for N_0 across runs, mean = ' num2str(mN0(1))]);
    legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    subplot(2, 1, 2);
    plot(1:M, relmse_r, 1:M, mr, 'k', 1:M, mr - 2*sr, 'r', 1:M, mr + 2*sr, 'r')
    xlabel('runs')
    ylabel('% (1 - rhat/r)^2')
    title(['Relative MSE(end) for r across runs, mean = ' num2str(mr(1))]);
    legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    
    
end