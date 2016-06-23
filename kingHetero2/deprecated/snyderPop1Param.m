% Perform Snyder filtering for for variable theta(t) = xsin(wt) + xmax

% Assumptions and Modifications
% - scaled time for the binfac idea
% - sinusoidal formulation for lam
% - the effective population varies with r parameters with m_i values
% - the data is k-coalescent times and the process self-exciting

clc
close all
clear all

%% Generate the appropriate coalescent data for a defined N(t)

% Set number of parameters and spaces
m = 100;
xmin = 0.1;
xmax = 10;
xset = linspace(xmin, xmax, m);
x = randsample(xset, 1);
w = 100;

% Set initial parameters: q0 is prior, N is true value, mi is dimension
% of parameter space, n is no. samples and m is total vector length
q0 = ones(1, m)/m;
n = 1000;
nData = n-1;

% Get coalescent binomial factors for lineages and set function type for
% rate ----- not the same as the function for N
nset = n:-1:2;
fac = nset.*(nset-1)/2;
%fac = ones(size(fac));


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
[twait, tcoal] = getCoalData1Param(fac, x, nData, w, xmax, 1); 
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
options = odeset('NonNegative', 1:m);
elemLen = zeros(1, nData);
xev = zeros(1, nData);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain approproate binomial factor
    binfac = fac(i);
    
    % Solve linear ODEs continuously with setting of options, no Q
    [tsol, qsol] = ode113(@(ts, y) odeSnyLinNonHomo1Param(ts, y, [], binfac, xset, w, 1),...
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
        % Sine function of time
        lamtdiag(j, 1:m) = binfac*(xset*sin(w*tsol(j)) + xmax);
    end
    
    % Perturb the q posterior for the new event
    lampert = diag(lamtdiag(end, :));
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :)*lampert./sum(qev(i+1, :)*lampert);
    
    % Estimate the parameters of the rate
    xhatset{i} = zeros(1, length(tsol));
    lamset{i} = zeros(1, length(tsol));
    thetaset{i} = zeros(1, length(tsol));
    for j = 1:length(tsol)
        % Parameters, rate and 1/N(t)
        lamset{i}(j) = qsol(j, :)*lamtdiag(j, :)';
        thetaset{i}(j) = (1/binfac)*qsol(j, :)*lamtdiag(j, :)';      
        xhatset{i}(j) = qsol(j, :)*xset';
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
end

%% Convert data in MSE estimates and plot comparisons

% Calculate true non-lineage dependent rate 
theta = x*sin(w*tn) + xmax;
tnmin = min(tn);
tnmax = max(tn);

lam = zeros(size(tn));
% Calculate poisson rate for coalescent events %% <<<<<< check for repeats
% in values due to duplicate tn values
for i = 1:nData
    id = find(tn < tcoal(i+1) & tn >= tcoal(i));
    tni = tn(id);
    lam(id) = fac(i)*(x*sin(w*tni) + xmax);
end

% Get the maximum and mean rate - by entering 0 one gets lam stats
[lamstat, lam_m, tm] = coalMSE(tn, lam, zeros(size(lam)));
lam_mean = lamstat.val(1);
lam_max = max(lam_m);
r = lam_max/lam_mean;
disp(['Log of max to mean rate is ' num2str(log(r))]);

% Plot the parameter and theta = 1/N estimates and rate
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

% % MSE calculation and time scales
% [statSet, em, tm] = coalMSE(tn, N, Nhat);
% tnmin = min(tn);
% tnmax = max(tn);
% 
% % Plot the estimate against N
% figure;
% plot(tn, Nhat, tn, N, 'k');
% xlabel('time');
% ylabel('effective population');
% legend('Nsny', 'N', 'location', 'best');
% title(['Snyder estimate for [n m] = [' num2str(n) ' ' num2str(m) ']']);
% 
% % Get squared error ratios
% RsqErr = (1 - Nhat./N).^2;
