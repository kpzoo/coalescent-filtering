%% Perform Snyder filtering for multi-variable N(t) on real data - no N

% Assumptions and Modifications
% - includes posterior for N(t)
% - allows simulation of log of Pybuse piecewise exponential
% - starts to use MBE toolbox
% - added confidence bounds and true Nt with last posterior, removed others
% - NEED TO ACCOUNT FOR DUPLICATE TIMES
% - calculation of tcoal using loops, hopefully a correction
% - removed plots including N and lam as no such values with real data
% - modified from batchMV2 to use real fasta sequence data
% - removed sinusoidal uniform sampling code
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
%clear all

% Booleans to control functionality and tree 
tic;
plotindiv = 1;
getTree = 0;
linBool = 0;

% Variables to control MBE use and dist type and duplicate treatment
useMBE = 0;
remDuplicates = 1;
dtype = 0;
readPybusxls = 0;

% Use the MBE toolbox
if useMBE
    % Add the path of the MBE toolbox
    cd ..
    ss = cd;
    ss = [ss '/MBEToolbox'];
    cd('HCV analysis');
    addpath(genpath(ss));
end
    
% Functional forms
fnid = 6;
global fnstr;
fnstr = {'x_1e^-x_2t', 'x_1sin(x_2t + x_3) + x_4', 'x_1',...
    'logistic growth', 'piecewise exponential', 'piecewise exponential', 'log piecewise'};
fnstr4 = 'x_1(1 + e^(-x_2x_3))/(1 + e^(-x_2(x_3 - t)) + x_4';
fnstr5 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_4) + x_1e^(-x_2(x_4 - x_3))I(x_4)';
fnstr6 = 'x_1I(x_3) + x_1e^(-x_2(t - x_3)I(x_3, x_3+x_4) + x_1e^(-x_2x_4)I(x_3+x_4)';
fnstr7 = 'x_1I(x_3) + (x_1 -x_2(t - x_3))I(x_3, x_3+x_4) + (x_1 -x_2x_4)I(x_3+x_4)';
numRVset = [2 4 1 4 4 4 4];

%% Extract n and the coalescent times from real fasta data

% Mutation rate clock and data from Pybus 2003
mu = 0.79*(10^-3);
name = 'result68.fasta';

% Main code to draw out coalescent data
if ~readPybusxls
    if ~useMBE
        % Coalescent times from UPGMA tree and Jukes-Cantor
        [tcoal, twait, trees, n, evDist, dcomp] = getCoalFasta(mu, name);
        if getTree
            % Sequence tree
            phytreeviewer(trees.treeSeq);
            % Coalescent time tree
            phytreeviewer(trees.treeCoal);
        end
    else
        % Coalescent times from UPGMA and MBE distance measures
        %%% <--------------- issue with length(tcoal) derived and zeros
        dbranch = getMBEdata(name, dtype);
        n = length(dbranch);
        tcoal = dbranch/mu;
        tcoal = tcoal';
        twait = diff(tcoal);
    end
else
    % Get coalescent times directly from xls file from Genie
    if strcmp(name, 'result63.fasta')
        coalData = xlsread('coalPybus63');
    else
        coalData = xlsread('coalPybus68');
    end
    twait = coalData(:, 2);
    tcoal = cumsum([0; twait]);
    n = length(tcoal);
end

% Other way to account for coalescent duplicates instead of maintaining
% the filter coalescent
if remDuplicates
    if(any(diff(tcoal) == 0))
        tcoal = unique(tcoal);
        nold = n;
        n = length(tcoal);
        disp(['Duplicate coalescents means n goes from ' num2str(nold) ' to ' num2str(n)]);
    end
end

%% Parameter settings

% Set number of parameters, dimensions (mi) and data length
mi = [20 20 20 20];
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
        minSpace = [1000 0 50 50];
        maxSpace = [10000 10 max(tcoal) 500];
    case 5
        % Pybus piecewise exponential
        minSpace = [5000 0 0 0];
        maxSpace = [15000 1 200 200];
    case 6
        % Pybus piecewise exponential but specify x and y-x vs y
        minSpace = [5000 0 30 0];
        maxSpace = [30000 0.75 120 80];
    case 7
        % Pybuse piecewise exponential but in terms of log(NC)
        minSpace = [log(5000) 0 0 0];
        maxSpace = [log(30000) 0.75 80 80];
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

% Get the random variable sets
xset = cell(1, 1);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
end
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

%% Filter the coalescent times from the fasta data

% Main function to perform Snyder estimation and return MMSE estimates
[xhat, Nhat, lamhat, tn, nLin, statsStruc, qn] = snyderRealData(fn, q0, tcoal);
tnmin = min(tn);
tnmax = max(tn);
tnlen = length(tn);
paramhat = xhat(end, :)';
fn.paramhat = paramhat;

% Extract outputs
xhatmean = statsStruc.xhatmean;
xstdlb = statsStruc.xstdlb;
xstdub = statsStruc.xstdub;

% Calculate a smooth function based on final value of parameters - pretend
% these are true values - this is suboptimal of course <---------------
t = linspace(tcoal(1), tcoal(end), 1000);
[Nt, ~] = getTimeVaryingN(fn, t, 1, paramhat);

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

% Reproduce Pybus2003 plot with time reorder and confidence intervals
if fnid == 5 || 6
    % Get other parameters of interest
    pm = [paramhat; 0];
    if fnid == 6
        % Param 4 is y-x
        pm(5) = paramhat(1)*exp(-paramhat(2)*(paramhat(4)));
    else
        % Param 4 is y
        pm(5) = paramhat(1)*exp(-paramhat(2)*(paramhat(4) - paramhat(3)));
    end
    infR = pm(1)/pm(5);
    
    % Get time in years and reorder so towards the present 
    tpres = 1993;
    maxPast = tcoal(end);
    tz = linspace(tpres-maxPast, tpres, 1000);
    tpast = 1895;
    
    % Parameter estimates from Pybus2003
    if strcmp(name, 'result68.fasta')
        % Data from 68 sequences
        oli.NC = [4095 10310 18960];
        oli.r = [0.075 0.264 0.620];
        oli.x = [1941 1953 1966];
        oli.y = [1924 1934 1943];
        oli.NA = [153 245 345];
        oli.TMRCA = [1258 1374 1481];
    else
        % Data from 63 sequences
        oli.NC = [3323 8779 15780];
        oli.r = [0.072 0.237 0.564];
        oli.x = [1941 1953 1966];
        oli.y = [1922 1932 1940];
        oli.NA = [99.6 170 251];
        oli.TMRCA = [1673 1710 1747];
    end
    
    % Get population optimal and bounds from last posterior
    idq = round(linspace(1, length(tz), 20));
    [Nf, Nbnd] = getPopForward(tz, qnlast, fn, tpres);
    %[Nf, Nbnd, qs, Nqs] = getPopForward2(tz, qnlast, fn, tpres, idq);
    
    % Get TMRCA from tree back in years
    TMRCA.yr = round(tpres - maxPast);
    TMRCA.len = maxPast;
    disp(['TMRCA = ' num2str(TMRCA.yr)]);
    
    % PAT treatment period
    patx = [1920 1920 1980 1980];
    Nmin = min(min(Nbnd));
    Nmax = max(max(Nbnd));
    paty = [Nmin Nmax Nmax Nmin];
    
    % Get the individual parameters in this form from marginal distribs
    est.NC = [xstdlb(end, 1), xhatmean(end, 1), xstdub(end, 1)];
    if fnid == 7
        est.NC = exp(est.NC);
    end
    est.NC = round(est.NC);
    est.r = [xstdlb(end, 2), xhatmean(end, 2), xstdub(end, 2)];
    est.x = tpres - [xstdub(end, 3), xhatmean(end, 3), xstdlb(end, 3)];
    est.x = round(est.x);
    if fnid == 5
        est.y = tpres - [xstdub(end, 4), xhatmean(end, 4), xstdlb(end, 4)];  
    else
        % Used deviation on y as deviation on y-x + mean x (approx)
        est.y = xhatmean(end, 3) + [xstdub(end, 4), xhatmean(end, 4), xstdlb(end, 4)];
        est.y = tpres - est.y;
    end
    est.y = round(est.y);
    
    % NA estimates more complicated - max when it is xmin - ymax and rmax
    % but NCmin with x, y defined with tpres
    if fnid == 5
        NAmean = est.NC(2)*exp(-est.r(2)*(est.x(2) - est.y(2)));
        xy = [est.x(1) - est.y(3) est.x(3) - est.y(1)];
        NAmax = est.NC(3)*exp(-est.r(1)*(min(xy)));
        NAmin = est.NC(1)*exp(-est.r(3)*(max(xy)));
    else
        NAmean = est.NC(2)*exp(-est.r(2)*(xhatmean(end, 4)));
        xy = [est.x(1) - est.y(3) est.x(3) - est.y(1)];
        NAmax = est.NC(3)*exp(-est.r(1)*(min(xy)));
        NAmin = est.NC(1)*exp(-est.r(3)*(max(xy)));
    end
    est.NA = [NAmin NAmean NAmax];
    est.NA = round(est.NA);
    clear NAmax NAmean NAmin
    
    % Convert coalescent times to years and get a vector to overlay on PAT
    if size(tcoal, 2) == 1
        % Set tcoal to a row vector if it isn't already
        tcoal = tcoal';
    end
    yr_coal = tpres - tcoal;
    n_coal = ones(1, n);
    vx_coal = repmat(yr_coal, 2, 1);
    v_coal = [zeros(1, n); xstdub(end, 1)*ones(1, n)];
    
    % Get space for NA from other spaces
    rmin = minSpace(2);
    rmax = maxSpace(2);
    xmin = minSpace(3);
    xmax = maxSpace(3);
    ymin = minSpace(4);
    ymax = maxSpace(4);
    NCmin = minSpace(1);
    NCmax = maxSpace(1);
    NAmax = NCmax*exp(-rmin*(ymin - xmax));
    NAmin = NCmin*exp(-rmax*(ymax - xmin));
end

%% Visualise data - no MSE calculation possible

% Plot best suboptimal N (assuming last value is best parameter value)
figure;
semilogy(t, Nt, 'k', 'LineWidth', 2);
hold on
semilogy(tn, Nhat, 'ro');
hold off
legend('final subst', 'iterative', 'location', 'best');
xlabel('time');
ylabel(['N(t) = ' fn.name]);
title(['N(t) final estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
xlim([tnmin tnmax]);
    
% Plot the N estimates and coalescent rate lambda
if plotindiv
    % The coalescent rate which includes lineages
%     figure;
%     semilogy(tn, lamhat);
%     xlabel('time');
%     ylabel('\lambda(x, t)');
%     title(['Rate estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
%     xlim([tnmin tnmax]);
    
%     % The population function
%     figure;
%     semilogy(tn, Nhat);
%     xlabel('time');
%     ylabel(['N(t) = ' fn.name]);
%     title(['N(t) estimate for [n m #params] = [' num2str(n) ' ' num2str(m) ' ' num2str(numRV) ']']);
%     xlim([tnmin tnmax]);
    
    % The parameters estimated and the std deviation bounds
    figure;
    for i = 1:numRV
        subplot(numRV, 1, i);
        plot(tn, xhatmean(:, i), 'k', tn, xstdlb(:, i),...
            'b', tn, xstdub(:, i), 'b');
        xlabel('time');
        ylabel(['x_' num2str(i) ' in N(t) =' fn.name]);
        title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
        legend('estimate', '-2 std', '+2 std', 'location', 'best');
        xlim([tnmin tnmax]);
    end
    
%     % Plot number of lineages across time
%     figure;
%     stairs(tcoal, n:-1:1);
%     xlabel('time');
%     ylabel('no. of lineages');
%     title('Number of lineages with time');
%     xlim([tnmin tnmax]);
    
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
        title(['Evolution of posteriors: param est = ' num2str(paramhat(i))]);
    end
    
    % Pybus2003 plot type
    if fnid == 5 || 6
        figure;
        fill(patx, paty, 'c:');
        hold on
        plot(tz, Nf, 'k', 'LineWidth', 2);
        plot(tz, Nbnd, 'r', 'LineWidth', 2);
        plot(vx_coal, v_coal, 'k:');
        hold off
        legend('PAT period', 'N(t)', 'lower bound', 'upper bound', 'location', 'best');
        xlim([tpast-100 tpres]);
        xlabel('year');
        ylabel('estimated effective number of infections');
        title('Estimated demographics of Egyptian HCV');
        
        % Another visualisation of the bound in relative terms
        NN = Nbnd./repmat(Nf, [2 1]);
        figure;
        plot(tz, NN);
        xlim([tpast-100 tpres]);
        xlabel('year');
        ylabel('bound(N(t))/N(t)');
        title('Relative confidence in N(t) for Egyptian HCV');
    end
end

% Clock time set to minutes
tIter = toc;
tIter = tIter/60;
disp(['Simulation time is ' num2str(tIter(end))]);
disp('********************************************************************');
for i = 1:numRV
    disp(['Estimate x_' num2str(i) ' = ' num2str(paramhat(i))]);
end
if fnid == 5 || 6
    % Additional dependent parameter calculated and ratio of infections
    disp(['Estimate x_5 = ' num2str(pm(5))]);
    disp(['Infection ratio = ' num2str(infR)]);
    if fnid == 5
        disp(['Length of exponential growth = ' num2str(pm(4) - pm(3))]);
    else
        disp(['Length of exponential growth = ' num2str(pm(4))]);
    end
end
    
% Store individual run removing heavy, unneeded variables
save(['real_' num2str(numRV)], 'fn', 'nLin', 'tn', 'lamhat', 'Nhat', 'statsStruc',...
    'tIter', 'tcoal', 'useMBE', 'remDuplicates', 'readPybusxls');