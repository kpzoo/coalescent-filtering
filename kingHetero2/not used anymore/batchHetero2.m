%% Perform Snyder filtering for multi-variable N(t)

% Assumptions and Modifications
% - removed tree option as not developed for heterochronous data
% - set for heterochronous data
% - includes fnstr6 which is reparametrised <-------- TO DO
% - works with getPopBackward to get optimal Nt estimates
% - includes N(t) from the final estimate of the parameters and bounds
% - the tn time vector is not uniformly sampled and has repetitions
% - uses x notation and expects certain functional forms
% - added repeatability for M runs
% - used rejection sampling to simulate NHPP
% - the data is k-coalescent times and the process self-exciting

clc
close all
clearvars

% Booleans to control functionality and runs, M and optimal N(t)
tic;
M = 1;
plotindiv = 1;
profiling = 0;
optimalNt = 1;

if profiling
    profile on
else
    profile off
end

% Time logging, boolean for MSE calc
tIter = zeros(1, M);
solveMSE_t = 0;

% Boolean structure for snyder ODE solution - linearity and single lineages
bools.linBool = 1;
bools.singLinRem = 1;

% Functional forms
fnid = 1;
global fnstr;
fnstr = {'x_1e^-x_2t', 'x_1sin(x_2t + x_3) + x_4', 'x_1',...
    'logistic', 'pybus exp'};
fnstr4 = 'x_1(1 + e^(-x_2x_3))/(1 + e^(-x_2(x_3 - t)) + x_4';
fnstr5 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_4) + x_1e^(-x_2(x_4 - x_3)I(x_4)';
numRVset = [2 4 1 4 4];

% Set heterochronous samples and lineages
svec = [0 0.5 1 1.5 2]*10^4;
svec = 0;
nvec = [20 20 20 20 20];
nvec = 200;
n = sum(nvec);
nSampTimes = length(svec);

% Set number of parameters, dimensions (mi) and data length
mi = 10*ones(1, 4);
if fnid == 1
    mi = [20 20];
end
numRV = length(mi);
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
fn.numRV = numRV;
fn.nEstRV = nEstRV;
fn.nSampTimes = nSampTimes;

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

% Generate thinning variable based on function
facSort = sort([0 fac]);
szf = size(facSort);
switch(fnid)
    case 1
        % Exponential
        Lmax = facSort;
    case 2
        % Sinusoidal
        Lmax = (fac(1)/(x(4) - x(1)))*ones(szf);
    case 3
        % Constant
        Lmax = fac(1)*ones(szf);
    case 4
        % Logistic functions
        Lmax = facSort/x(4);
    case 5
        % Pybus piecewise exponential - x3 and x4 are both times 
        Lmax = facSort/(x(1)*exp(-x(2)*(x(4) - x(3))));
end

% Extra function entries for heterochronous sampling
fn.svec = svec;
fn.nvec = nvec;
fn.Lmax = Lmax;
fn.facSort = facSort;

%% Main Loop

for ii = 1:M
    
    %% Generate coalescent data for a defined N(t)
    
    % Generate the coalescent data using rejection sampling
    events = simHeteroRejSampling(fn);
    tLin = events.tLin;
    nLin = events.nLin;
    
    % Lineage through time plot with sample times
    if plotindiv
        figure;
        stairs(tLin, nLin);
        hold on
        plot(repmat(svec, [2 1]), [zeros(size(svec)); n*ones(size(svec))], 'r--', 'linewidth', 2);
        hold off
        legend('lineages', 'sample times', 'location', 'best');
        xlabel('time');
        ylabel('no. lineages');
        title(['Lineages through time with no. sample times = ' num2str(nSampTimes)]);
        xlim([tLin(1) tLin(end)]);
        ylim([min(nLin) max(nLin)]);
    end
    
    % Get waiting times and tree constructions
    tcoal = events.tcoal;
    twait = diff(tcoal);
    
    %% Filter the simulated event times from the appropriate Poisson process
    
    % Main function to perform Snyder estimation and return MMSE estimates
    [xhat, Nhat, lamhat, lam, N, tn, statsStruc, qn] = snyderHetero2(fn, q0, events, bools);
    tnmin = min(tn);
    tnmax = max(tn);
    tnlen = length(tn);
    
    % Marginalise the posterior across sim times - take 20 of them
    nPoster = 20;
    idpost = round(linspace(1, length(tn), nPoster));
    qmarg = cell(1, nPoster);
    probSums = cell(1, nPoster);
    for i = 1:nPoster
        [qmarg{i}, probSums{i}] = marginalise(numRV, IDMx, qn(idpost(i), :), mi);
    end
    qnlast = qn(end, :);
    clear qn
    
    % Reorganise into probability matrices for each variable
    qnp = cell(1, numRV);
    for i = 1:numRV
        qnp{i} = zeros(nPoster, mi(i));
        for j = 1:nPoster
            qnp{i}(j, :) = qmarg{j}{i};
        end
    end
    
    % Extract outputs
    xhatmean = statsStruc.xhatmean;
    xstdlb = statsStruc.xstdlb;
    xstdub = statsStruc.xstdub;
    
    % Display structure for parameters
    est.x = x;
    est.xhat = xhat(end, :);
    
    % Get optimal estimated population and bounds from last posterior
    if optimalNt
        tz = linspace(tcoal(1), tcoal(end), 1000);
        % Final estimates - only works for inhomogeneous MCs
        [Nf, Nbnd] = getPopBackward(tz, qnlast, fn);
        % True population
        [Nft, ~] = getTimeVaryingN(fn, tz, -1, x');
    end
    
    %% Convert data in MSE estimates and plot comparisons
    
    % MSE ratio calculation as a percent across time
    if solveMSE_t
        [Nstat, ~, ~, setData] = coalMSE(tn, Nhat./N, ones(size(tn)));
        R = 100*Nstat.val(3);
    else
        R = -1;
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
        
        % Posteriors and final estimates
        figure;
        for i = 1:numRV
            subplot(numRV, 1, i);
            plot(xset{i}, qnp{i}(1, :), 'b-', xset{i}, qnp{i}(end, :), 'ko-',...
                xset{i}, qnp{i}(2:end-1, :), 'g--');
            hold on
            plot([paramhat(i) paramhat(i)], [0 max(max(qnp{i}))], 'k');
            hold off
            % Condition to account for variables with only one value
            if min(xset{i}) ~= max(xset{i})
                xlim([min(xset{i}) max(xset{i})]);
            end
            xlabel(['x_' num2str(i)]);
            ylabel(['P(x_' num2str(i) '|data)']);
            legend('prior', 'posterior', 'intermediates', 'location', 'best');
            title(['Evolution of posteriors: param est = ' num2str(xhat(end, i))]);
        end
        
        
        % The parameters estimated and the std deviation bounds
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
        
        % Joint final posterior if bivariate
        if numRV == 2 && nEstRV == 2
            % Reshape last posterior to 2D form for bivariate plotting
            qtemp = reshape(qnlast, mi);
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
        
        % Final N(t) estimate from final parameters with coalescent points
        if optimalNt
            vx_coal = repmat(tcoal, 2, 1);
            v_coal = [zeros(1, n); max(max(Nbnd))*ones(1, n)];
            vx_samp = repmat(svec(2:end), 2, 1);
            v_samp = [zeros(1, nSampTimes-1); max(max(Nbnd))*ones(1, nSampTimes-1)];
            figure;
            plot(tz, Nft, 'r', 'LineWidth', 2);
            hold on
            plot(tz, Nf, 'k', 'LineWidth', 2);
            plot(tz, Nbnd, 'g--', 'LineWidth', 2);
            plot(vx_samp, v_samp, 'm');
            plot(vx_coal, v_coal, 'k:');
            hold off
            legend('N(t)', 'Nhat(t)', 'lower bound', 'upper bound', 'sample times', 'coalescents', 'location', 'best');
            xlim([tcoal(1) tcoal(end)]);
            xlabel('time');
            ylabel('N(t) estimate and bounds');
            title('Population estimate using final posterior');
            grid;
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
    save(['batchHet_' num2str(numRV) '_' num2str(M)], 'mse_t', 'param', 'paramhat',...
        'n', 'M', 'tIter', 'fn', 'svec', 'nvec');
    
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
        grid;
    end

else
    % Store individual run removing heavy, unneeded variables
    save(['hetero_' num2str(numRV)], 'fn', 'tz', 'Nf', 'Nft', 'Nbnd', 'N',...
        'Nhat', 'statsStruc', 'mse_t', 'tIter', 'est', 'events');
end