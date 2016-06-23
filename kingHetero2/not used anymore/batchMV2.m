%% Perform Snyder filtering for multi-variable N(t)

% Assumptions and Modifications
% - modified from batchMV to include piecewise exp
% - removed minSpace and maxSpace as cells of all fns to single fn arrays
% - the tn time vector is not uniformly sampled and has repetitions
% - added code to calculate a matrix of identifiers for function calling
% - removed product calculation on parameter estimates
% - uses x notation and expects certain functional forms
% - added repeatability for M runs
% - assumes appropriate data provided
% - used rejection sampling to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clear all

% Booleans to control functionality and runs, M
tic;
M = 1;
plotindiv = 1;
profiling = 0;


if profiling
    profile on
else
    profile off
end

% Time logging variable and boolean to control MSE calc and draw tree
tIter = zeros(1, M);
solveMSE_t = 0;
getTree = 0;

% Functional forms
fnid = 2;
global fnstr;
fnstr = {'x_1e^-x_2t', 'x_1sin(x_2t + x_3) + x_4', 'x_1',...
    'logistic growth', 'piecewise exponential'};
fnstr4 = 'x_1(1 + e^(-x_2x_3))/(1 + e^(-x_2(x_3 - t)) + x_4';
fnstr5 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_4) + x_1e^(-x_2(x_4 - x_3)I(x_4)';
numRVset = [2 4 1 4 4];

% Set number of parameters, dimensions (mi) ad data length
mi = [8 8 8 8]*2.5/2;
numRV = length(mi);
n = 1000;
nData = n-1;
m = prod(mi);

% Check for pre-specified parameters not to be estimated - mi values of 1
specParam = find(mi == 1);
nEstRV = numRV - length(specParam);

% Get coalescent binomial factors for lineages and change notation
nset = n:-1:2;
fac = nset.*(nset-1)/2;

% Space sets for functions
switch(fnid)
    case 1
        % Exponential
        minSpace = [100*n 0.1];
        maxSpace = [1000*n 10];
    case 2
        % Sinusoidal
        minSpace = [100*n 0.1 0 1100*n];
        maxSpace = [1000*n 10 pi/2 1200*n];
    case 3
        % Constant
        minSpace = 100*n;
        maxSpace = 1000*n;
    case 4
        % Logistic
        minSpace = [100*n 0.1 1000 50*n];
        maxSpace = [1000*n 10 5000 100*n];
    case 5
        % Pybus piecewise exponential
        minSpace = [1000*n 0.1 5 40];
        maxSpace = [2000*n 0.75 30 70];
end

% Modify spaces if specify parameters so that they are set to min values
if nEstRV > 0
    minSpace(specParam) = maxSpace(specParam);
end

% Function identifier structure
fn.id = fnid;
fn.name = fnstr(fnid);
fn.nData = nData;
fn.mi = mi;
fn.m = m;
fn.fac = fac;
fn.numRV = numRV;
fn.nEstRV = nEstRV;

%% Initialise and define key parameters and storage variables

% Check numRV compatible with fntype
if numRV ~= numRVset(fnid)
    error('Initial settings on RV space inconsistent');
end

% Set uniform prior q0
q0 = ones(1, m)/m;

% Get the instance of the random variables for rate function, lam
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    x(i) = datasample(xset{i}, 1);
end

% Parameters to store results with MSE across time, J_t = mse_t
mse_t = zeros(1, M);
param = zeros(numRV, M);
paramhat = zeros(numRV, M);
fn.param = x;
fn.xset = xset;

% Create a matrix of identifiers to tell which xset{i} values are used for
% each entry in N(t) and lam(t) calculations
IDMx = zeros(numRV, m);
% Initialise with first variable which has no element repetitions
idxset = 1:mi(1);
IDMx(1, :) = repmat(idxset, 1, m/mi(1));
for i = 2:numRV
    % For further variables numReps gives the number of set repetitions
    % while kronVec gives the number of element repetitions
    idxset = 1:mi(i);
    numReps = m/prod(mi(1:i));
    kronVec = ones(1, prod(mi(1:i-1)));
    IDMx(i, :) = repmat(kron(idxset, kronVec), 1, numReps);
end

% Get the values corresponding to the matrix
xsetMx = zeros(numRV, m);
for i = 1:numRV
    xsetMx(i, :) = xset{i}(IDMx(i, :));
end
fn.IDMx = IDMx;
fn.xsetMx = xsetMx;

for ii = 1:M
    
    %% Generate coalescent data and tree for a defined N(t)
    
    % Generate the coalescent data using rejection sampling
    [twait, tcoal] = getCoalDataFn(fn);
    C = cumsum(twait);
    
    % Generate coalescent tree
    if getTree
        % Coalescent time based tree
        [children, labels, tree] = getKingmanTree(n, C');
        % Gives sparse connection matrix of tree
        [M, ID, Dists] = getmatrix(tree);
        phytreeviewer(tree);
    end
    
    %% Filter the simulated event times from the appropriate Poisson process
    
    % Main function to perform Snyder estimation and return MMSE estimates
    [xhat, Nhat, lamhat, lam, N, tn, nLin, statsStruc] = snyder(fn, q0, tcoal);
    tnmin = min(tn);
    tnmax = max(tn);
    tnlen = length(tn);
    
    % Extract outputs
    xhatmean = statsStruc.xhatmean;
    xstdlb = statsStruc.xstdlb;
    xstdub = statsStruc.xstdub;
    
    % Display structure for parameters
    est.x = x;
    est.xhat = xhat(end, :);

    % Twenty iterations of the posterior from prior
    if numRV == 2
        phat = statsStruc.phat;
        qtemp = statsStruc.qtemp;
        nPoster = 20;
        qnp = cell(1, numRV);
        for kk = 1:numRV
            finidpost = tnlen;
            idpost = round(linspace(1, finidpost, nPoster));
            % Posteriors spaced as specified by idpost
            qnp{kk} = phat{kk}(idpost, :);
        end
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % MSE ratio calculation as a percent across time
    if solveMSE_t
        [Nstat, ~, ~, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
        R = 100*Nstat.val(3);
    else
        R = -1;
    end
    if fnid == 2 && plotindiv
        % Note - cannot resample for Nhat etc because var(tn) >> E[tn] due
        % to the lineages - would need to rescale somehow <----------
        
        % Resample data uniformly and use parameters if sin
        L = 10;
        F = 20;
        [tx, ~, ~] = makeUniformSamples(tn, L, F, Nhat, lamhat);
        Nx = x(1)*sin(x(2)*tx + x(3)) + x(4);
        Nhatx = xhat(end, 1)*sin(xhat(end, 2)*tx + xhat(end, 3)) + xhat(end, 4);
        % Plot results
        figure;
        plot(tx, Nx, tx, Nhatx);
        xlim([min(tx) tnmax]);
        xlabel('time');
        ylabel(['N(t) = ' fn.name]);
        legend('N', 'Nsubst', 'location', 'best');
        title(['N(t) suboptimal estimate for [n m #params] = [' num2str(n)...
            ' ' num2str(m) ' ' num2str(numRV) ']']);
    end
    
    % Plot the N estimates and coalescent rate lambda
    if plotindiv
        % The coalescent rate which includes lineages
        figure;
        semilogy(tn, lamhat, tn, lam, 'r');
        xlabel('time');
        ylabel('\lambda(x, t)');
        legend('lamsny', 'lam', 'location', 'best');
        title(['Rate estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
        % The population function
        figure;
        semilogy(tn, Nhat, tn, N, 'r');
        xlabel('time');
        ylabel(['N(t) = ' fn.name]);
        legend('Nsny', 'N', 'location', 'best');
        title(['N(t) estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
        xlim([tnmin tnmax]);
        
        % The parameters estimated and the std deviation bounds
        if numRV ~= 2
            figure;
            for i = 1:numRV
                subplot(numRV, 1, i);
                plot(tn, x(i)*ones(size(tn)), 'r', tn, xhatmean(:, i), 'k', tn, xstdlb(:, i),...
                    'b', tn, xstdub(:, i), 'b');
                xlabel('time');
                ylabel(['x_' num2str(i) ' in N(t) =' fn.name]);
                title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
                legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
                xlim([tnmin tnmax]);
            end
        else
            figure;
            for i = 1:numRV
                subplot(numRV, 1, i);
                plot(tn, x(i)*ones(size(tn)), 'r', tn, xhatmean{i}, 'k', tn, xstdlb{i},...
                    'b', tn, xstdub{i}, 'b');
                xlabel('time');
                ylabel(['x_1 in N(t) =' fn.name]);
                title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
                legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
                xlim([tnmin tnmax]);
            end
            
            % Posteriors and final estimates
            figure;
            for i = 1:numRV
                subplot(numRV, 1, i);
                plot(xset{i}, qnp{i}(1, :), 'b', xset{i}, qnp{i}(end, :), 'k',...
                    xset{i}, qnp{i}(2:end-1, :), 'g--');
                hold on
                plot([x(i) x(i)], [0 max(max(qnp{i}))], 'k');
                hold off
                xlabel(['x_' num2str(i)]);
                ylabel(['P(x_' num2str(i) '|data)']);
                legend('prior', 'posterior', 'intermediates', 'location', 'best');
                title(['Evolution of posteriors: [x_' num2str(i) ' x_' num2str(i) 'hat] = '...
                    [num2str(x(i)) ', ' num2str(xhat(end, i))]]);
            end
        end
        
        % Plot number of lineages across time
        figure;
        stairs(tn, nLin);
        xlabel('time');
        ylabel('no. of lineages');
        title('Number of lineages with time');
        xlim([tnmin tnmax]);
        
        % Joint final posterior if bivariate
        if numRV == 2 && nEstRV == 2
            if ~isnan(qtemp)
                figure;
                meshz(xset{1}, xset{2}, qtemp');
                xlabel('x_1');
                ylabel('x_2');
                zlabel('P(x_1, x_2 | data)');
                title(['Joint posterior with true [x1 x2] = ' num2str(x(1)) ' ' num2str(x(2))]);
            else
                warning('Mat:qNaN', 'qtemp has NaN values');
            end
        end
    end
    
    % Store batch results with param and paramhat having rows for numRVs
    mse_t(ii) = R;
    param(:, ii) = x';
    paramhat(:, ii) = xhat(end, :)';
    disp('***************************************************************');
    disp(['The MSE perc ratio is ' num2str(R)]);
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
    tIter(ii) = toc;
    
end

if profiling
    profile viewer
    profile off
end

% Clock time set to minutes
tIter = tIter/60;
disp(['Total time for ' num2str(M) ' iterations is ' num2str(tIter(end))]);
if M > 1
    tavgIter = mean(diff(tIter));
    disp(['Average time with [n m] = ' [num2str(n) ' ' num2str(m)] ' is ' num2str(tavgIter)]);

    % Save important bits of data
    save(['batch_' num2str(numRV) '_' num2str(M)], 'mse_t', 'param', 'paramhat',...
        'n', 'M', 'tIter', 'fn');
    
    % Parameter estimate errors in multivariate form
    eparam = param - paramhat;
    esq = eparam.^2;
    
    % Relative MSE across runs, J_r (from report)
    mse_r = 100*esq./(param.^2);
    mx = mean(mse_r, 2);
    sx = std(mse_r, 0, 2);
    
    % Plot the behaviour across runs
    figure;
    for i = 1:numRV
        subplot(numRV, 1, i);
        mxi = mx(i)*ones(1, M);
        sxi = mx(i)*ones(1, M);
        plot(1:M, mse_r(i, :), 1:M, mxi, 'k', 1:M, mxi - 2*sxi, 'r', 1:M, mxi + 2*sxi, 'r')
        xlabel('runs')
        ylabel('% (1 - xhat/x)^2')
        title(['Relative MSE(end) for x_' num2str(i) ' across runs, mean = ' num2str(mxi(1))]);
        legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    end

else
    % Store individual run removing heavy, unneeded variables
    save(['indiv_' num2str(numRV)], 'fn', 'nLin', 'tn', 'lamhat', 'lam', 'N',...
        'Nhat', 'statsStruc', 'mse_t', 'tIter', 'C', 'tcoal', 'est');
end