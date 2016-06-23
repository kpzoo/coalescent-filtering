%% Script that compares con-exp-con and exp functions sampled across time

% Assumptions and Modifications
% - removed neff and cap and lamz etc and many plots
% - compared both pybus and exponentials for same x1 and x2
% - removed plots, Nhatz, xhatvar etc and marginalisations
% - difference from script is input n, nSampTimes and maxSamp
% - only pybus piecewise exponential demographic functions
% - introduces the effective lineage count idea
% - in heterochronous case divides samples uniformly %% <----- can change to other distributions later

clc
close all
clearvars
tic; 

% Boolean to decide comparison
comparison = 1;

% Set number of parameters and size of discretised grid
numRV = 4;
mi = 20*ones(1, numRV);
if length(mi) ~= numRV
    error('mi is of incorrect length');
end
m = prod(mi);

% Set parameter space and true values but specify x and y-x vs y
%minSpace = [10000 0 50 20];
%maxSpace = [50000 0.1 75 50];

minSpace = [10000 0 50 50];
maxSpace = [100000 0.1 50 50];

% Get the instance of the random variables for demographic function
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    if comparison
        % Set to middle of range
        x(i) = xset{i}(mi(i)/2);
    else
        % Randomly sample
        x(i) = datasample(xset{i}, 1);
    end
end

% Set number of samples and assume uniformly divided in time
n = 200;
nSampTimes = 1;
nvec = round(n/nSampTimes)*ones(1, nSampTimes);
maxSamp = 2*x(3) + x(4);
if nSampTimes > 1
    svec = linspace(0, maxSamp, nSampTimes);
else
    % linspace picks maxSampTime if nSampTimes = 1
    svec = 0;
end

% Set no. data points as distinct coalescent events and sample times -
% becomes n-1 for isochronous when nSampTimes = 1
nData = n + nSampTimes - 2;

% Get coalescent binomial factors for lineages and make no. lineages an
% index so facSort(1) = 0
nset = n:-1:2;
fac = nset.*(nset-1)/2;
facSort = sort([0 fac]);

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

%% Generate the coalescent data with heterochronous sampling via rescaling

% Simulation of con-exp-con
[nLin, tLin, tcoal] = simExpPyb(svec, nvec, n, nData, x, 0, facSort, nSampTimes);
% Simulation of exponential
[nLin1, tLin1, tcoal1] = simExpPyb(svec, nvec, n, nData, x, 1, facSort, nSampTimes);

% Check for case when x5 goes to zero or tcoal is not distinct
if any(diff(tcoal) == 0) || any(diff(tcoal1) == 0)
    error('The coalescent times are not distinct');
end

% Get N(t) con-exp-con in time - true function
tz = linspace(0, 2*x(3) + x(4), 1000);
Nz = zeros(size(tz));
for i = 1:length(tz)
    I1 = tz(i) <= x(3);
    I2 = tz(i) > x(3) & tz(i) < x(4) + x(3);
    I3 = tz(i) >= x(4) + x(3);
    Nz(i) = x(1).*I1 + x(1).*exp(-x(2).*(tz(i) - x(3))).*I2 + ...
        x(1).*exp(-x(2).*x(4)).*I3;
end

% Get N(t) exponential in time - true function 
Nz1 = x(1).*exp(-x(2)*tz);


%% Snyder Filtering of the coalescent-heterochronous data

% Con-exp-con filtering
[xhat, xhatmean, xstdlb, xstdub, qnlast] = pybexpFilter(mi, nData, 0, xsetMx, tLin, nLin, facSort, svec);
% Display structure for parameters
est.x = x;
est.xhat = xhat(end, :);

% Final estimates for con-exp-con
fn.xsetMx = xsetMx;
fn.id = 6;
[Nzhat, Nbnd] = getPopBackward(tz, qnlast, fn);

% Exponential filtering
[xhat1, xhatmean1, xstdlb1, xstdub1, qnlast1] = pybexpFilter(mi, nData, 1, xsetMx, tLin1, nLin1, facSort, svec);
% Display structure for parameters
est1.x = x;
est1.xhat = xhat1(end, :);

% Final estimates for exponential
fn1.xsetMx = xsetMx(1:2, 1:prod(mi(1:2)));
fn1.id = 1;
[Nzhat1, Nbnd1] = getPopBackward(tz, qnlast1, fn1);

% Marginalise final posterior
[qmarg, probSums] = marginalise(numRV, IDMx, qnlast, mi);
[qmarg1, probSums1] = marginalise(2, IDMx(1:2, 1:prod(mi(1:2))), qnlast1, mi);

% Simulation time
tsim = toc/60;
disp(['Simulation time = ' num2str(tsim) ' mins']);

%% Plots to visualise comparative results

% Population functions
figure;
plot(tz, Nz, 'k', tz, Nzhat, 'r', tz, Nbnd, 'g--', 'linewidth', 2);
h = gca;
hold on
vx_samp = repmat(svec(2:end), 2, 1);
v_samp = [zeros(1, nSampTimes-1); h.YLim(2)*ones(1, nSampTimes-1)];
plot(vx_samp, v_samp, 'm');
plot(tz, Nz1, 'k', tz, Nzhat1, 'r', tz, Nbnd1, 'g--', 'linewidth', 2);
hold off
grid;
legend('N(t)', 'Nhat(t)', 'bounds', 'sample times', 'location', 'best');
xlim([tz(1) tz(end)]);
xlabel('time');
ylabel('population size');
title('Population estimates for con-exp-con and exponential');

% Compare posteriors for the first 
figure;
subplot(2, 1, 1);
hAx = plotyy(xset{1}, qmarg{1}, xset{1}, qmarg1{1});
hold all
plot([est.xhat(1) est.xhat(1)], hAx(1).YLim, 'b--');
plot([est1.xhat(1) est1.xhat(1)], hAx(2).YLim, 'r--');
hold off
legend('q(x_1) con-exp-con', '<x_1> con-exp-con', '<x_1> exponential', 'q(x_1) exponential', 'location', 'best');
xlabel('x_1 space');
grid;
title('Marginal posterior for x_1');
subplot(2, 1, 2);
hAx = plotyy(xset{2}, qmarg{2}, xset{2}, qmarg1{2});
hold all
plot([est.xhat(2) est.xhat(2)], hAx(1).YLim, 'b--');
plot([est1.xhat(2) est1.xhat(2)], hAx(2).YLim, 'r--');
hold off
legend('q(x_2) con-exp-con', '<x_2> con-exp-con', '<x_2> exponential', 'q(x_2) exponential', 'location', 'best');
xlabel('x_2 space');
grid;
title('Marginal posterior for x_2');
