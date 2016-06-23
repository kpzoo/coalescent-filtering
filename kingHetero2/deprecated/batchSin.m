% Perform Snyder filtering for variable N(t) = x1sin(x2)t + A with x1
% and x2 simultaneously estimated and A > x1max

% Assumptions and Modifications
% - similar to batchExp2Param1 but for sinusoid and keeping x notation
% - added repeatability for M runs
% - used rejection sampling to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clear all

% Booleans to control functionality
tic;
M = 1;
plotindiv = 1;

% Time logging variable and boolean to control MSE calc or not
tIter = zeros(1, M);
solveMSE_t = 0;

%% Initialise and define key parameters and storage variables

% Set number of parameters, dimensions (mi) abd data length
mi = [20 20];
numRV = length(mi);
n = 200;
nData = n-1;

% Set space, prior q0 and total dimension m
minspace = [100*n 0.01];
maxspace = [1000*n 0.1];
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

% Get (x1, x2) combined into a single array via kronecker products and
% define maximum offset
space = kron(xset{1}, xset{2});
A = 1.1*maxspace(1);

% Parameters to store results
MSEpercratio = zeros(1, M);
param = zeros(numRV, M);
paramhat = zeros(numRV, M);

for ii = 1:M
    
    %% Generate the appropriate coalescent data for a defined N(t)
    
    % Generate the coalescent data using rejection sampling
    [twait, tcoal] = getCoalDataSin(fac, x(1), nData, x(2), A);
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
    
    % Run Snyder filter across the time series
    for i = 1:nData
        % Obtain approproate binomial factor
        binfac = fac(i);
        
        % Solve linear ODEs continuously with setting of options, no Q
        [tsol, qsol] = ode113(@(ts, y) odeSnySin2Param(ts, y, [], binfac, xset, mi, A),...
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
            % Calculate lamt across time explicitly
            Nt = kron(xset{1}, sin(xset{2}*tsol(j))) + A;
            lamtdiag(j, :) = binfac./Nt;
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
    x1hat = phat{1}*xset{1}';
    x2hat = phat{2}*xset{2}';
    % Std deviation of conditional mean
    x1hatvar = phat{1}*(xset{1}'.^2) - x1hat.^2;
    x2hatvar = phat{2}*(xset{2}'.^2) - x2hat.^2;
    x1hatstd = sqrt(x1hatvar);
    x2hatstd = sqrt(x2hatvar);
    % Calculate upper and lower std deviation bounds mean +/- 2*std
    x1ubstd = x1hat + 2*x1hatstd;
    x1lbstd = x1hat - 2*x1hatstd;
    x2ubstd = x2hat + 2*x2hatstd;
    x2lbstd = x2hat - 2*x2hatstd;
    
    %% Convert data in MSE estimates and plot comparisons
   
    % Calculate true non-lineage dependent rate and population
    N = x(1)*sin(x(2)*tn) + A;
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
    
%     % Resample data uniformly and use parameters
%     L = 1000;
%     F = 20;
%     [tx, ~, ~] = makeUniformSamples(tn, L, F, Nhat, lamhat);
%     Nx = x(1)*sin(x(2)*tx) + A;
%     Nhatx = xhat(end, 1)*sin(xhat(end, 2)*tx) + A; %%<------ interpolation seems bad
    
    % Display structure
    xprodhat = prod(xhat, 2);
    xprod = prod(x);
    pmc = corr2(xhat(:,1), xhat(:,2));
    est.x1 = x(1);
    est.x1hat = xhat(end, 1);
    est.x2 = x(2);
    est.xhat2 = xhat(end, 2);
    est.xprod = xprod;
    est.xprodhat = xprodhat(end);
    est.pmc = pmc;
    
    % Plot the N estimates and coalescent rate lambda
    if plotindiv
        % The coalescent rate which includes lineages
        figure;
        semilogy(tn, lamhat, tn, lam, 'r');
        xlabel('time');
        ylabel('\lambda(x_1, x_2, t)');
        legend('lamsny', 'lam', 'location', 'best');
        title(['Rate estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
        % The population function
        figure;
        plot(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel('N(t) = x_1sin(x_2t) + A');
        legend('Nsny', 'N', 'location', 'best');
        title(['N(t) estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
%         % The parameters estimated
%         figure;
%         subplot(2, 1, 1);
%         plot(tn, x(1)*ones(size(tn)), 'r', tn, xhat(:, 1));
%         xlabel('time');
%         ylabel('x_1 in N(t) = x_1sin(x_2t) + A');
%         title(['Estimate of x_1 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
%         ylim([0.95*min([x(1); xhat(:, 1)]) 1.05*max([x(1); xhat(:, 1)])]);
%         xlim([tnmin tnmax]);
%         subplot(2, 1, 2);
%         plot(tn, x(2)*ones(size(tn)), 'r', tn, xhat(:, 2));
%         xlabel('time');
%         ylabel('x_2 in N(t) = x_1sin(x_2t) + A');
%         title(['Estimate of x_2 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
%         ylim([0.95*min([x(2); xhat(:, 2)]) 1.05*max([x(2); xhat(:, 2)])]);
%         xlim([tnmin tnmax]);
        
        % The parameters estimated and the std deviation bounds
        figure;
        subplot(2, 1, 1);
        plot(tn, x(1)*ones(size(tn)), 'r', tn, xhat(:, 1), 'k', tn, x1lbstd, 'b', tn, x1ubstd, 'b');
        xlabel('time');
        ylabel('x_1 in N(t) = x_1sin(x_2t) + A');
        title(['Estimate of x_1 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
        xlim([tnmin tnmax]);
        subplot(2, 1, 2);
        plot(tn, x(2)*ones(size(tn)), 'r', tn, xhat(:, 2), 'k', tn, x2lbstd, 'b', tn, x2ubstd, 'b');
        xlabel('time');
        ylabel('x_2 in N(t) = x_1sin(x_2t) + A');
        title(['Estimate of x_2 for [n m] = [' num2str(n) ' ' num2str(m) ']']);
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
        xlabel('x_1');
        ylabel('P(x_1|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: [x_1 x_1 hat] = ' [num2str(x(1)) ', ' num2str(xhat(end, 1))]]);
        subplot(2, 1, 2);
        plot(xset{2}, qnp{2}(1, :), 'b', xset{2}, qnp{2}(end, :), 'k',...
            xset{2}, qnp{2}(2:end-1, :), 'g--');
        hold on
        plot([x(2) x(2)], [0 max(max(qnp{2}))], 'k');
        hold off
        xlabel('x_2');
        ylabel('P(x_2|data)');
        legend('prior', 'posterior', 'intermediates', 'location', 'best');
        title(['Evolution of posteriors: [x_2 x_2 hat] = ' [num2str(x(2)) ', ' num2str(xhat(end, 2))]]);
        
        % Joint final posterior
        figure;
        meshz(xset{1}, xset{2}, qtemp);
        xlabel('x_1');
        ylabel('x_2');
        zlabel('P(x_1, x_2 | data)');
        title(['Joint posterior with true [x1 x2] = ' num2str(x(1)) ' ' num2str(x(2))]);
        
        % Product of parameters estimation
        figure;
        plot(tn, xprod*ones(size(tn)), tn, xprodhat);
        xlabel('time');
        ylabel('x_1*x_2 in N(t) = x_1sin(x_2t) + A');
        title(['Estimate of x_1*x_2 for [n m pmcc] = [' num2str(n) ' ' num2str(m) ' ' num2str(pmc) ']']);
        legend('true', 'estimate', 'location', 'best');
        xlim([tnmin tnmax]);
    end
    
    % MSE ratio calculation as a percent
    if solveMSE_t
        [Nstat, ~, ~, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
        R = 100*Nstat.val(3);
    else
        R = -1;
    end
    
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
    relmse_x1 = 100*esq(1, :)/(x(1)^2);
    relmse_x2 = 100*esq(2, :)/(x(2)^2);
    mx1 = mean(relmse_x1)*ones(1, M);
    sx1 = std(relmse_x1)*ones(1, M);
    mx2 = mean(relmse_x2)*ones(1, M);
    sx2 = std(relmse_x2)*ones(1, M);
    
    % Plot the behaviour across runs
    figure;
    subplot(2, 1, 1);
    plot(1:M, relmse_x1, 1:M, mx1, 'k', 1:M, mx1 - 2*sx1, 'r', 1:M, mx1 + 2*sx1, 'r')
    xlabel('runs')
    ylabel('% (1 - Nhat/N_0)^2')
    title(['Relative MSE(end) for x_1 across runs, mean = ' num2str(mx1(1))]);
    legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    subplot(2, 1, 2);
    plot(1:M, relmse_x2, 1:M, mx2, 'k', 1:M, mx2 - 2*sx2, 'r', 1:M, mx2 + 2*sx2, 'r')
    xlabel('runs')
    ylabel('% (1 - rhat/r)^2')
    title(['Relative MSE(end) for x_2 across runs, mean = ' num2str(mx2(1))]);
    legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    
    
else
    % Store individual run removing heavy, unneeded variables
    varlist = {'qset', 'qsol', 'qtemp', 'options', 'lamtdiag', 'lamnoTdiag',...
        'qn', 'qev', 'qnoev', 'lampert'};
    clear(varlist{:});
    save('tempIndiv');
end