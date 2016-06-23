% Function to perform Snyder filtering and estimation
function [xhat, Nhat, lamhat, lam, N, tn, statsStruc, qn] = snyderHetero(fn, q0, events, bools)

% Assumptions and modifications
% - nData is no longer n-1 but includes sample time events
% - modified for heterochronous sampling
% - removed nLin calculation from this function
% - allows linear vs non-linear
% - outputs qn and removes special numRV == 2 case
% - works for multivariable functions defined in fn
% - includes new functions for calculating lamdiag and N
% - uses matrices xsetMx etc defined in fn

% Extract inputs about fn
facSort = fn.facSort;
m = fn.m;
xsetMx = fn.xsetMx;
x = fn.param;
linBool = bools.linBool;

% Extract heterochronous inputs
tLin = events.tLin;
nLin = events.nLin;
singLinRem = bools.singLinRem;

% Data here is coalescent and sample events % <--------------------
nData = length(tLin)-1;

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
%options = odeset('NonNegative', 1:m);
options= odeset('RelTol', 1e-9, 'NonNegative', 1:m);
elemLen = zeros(1, nData);
lamtdiag = cell(1, 1);

% Loop across no. coalescents and no. sample times
nLoop = length(nLin);

% Run Snyder filter across the time series
for i = 1:nLoop-1
    % Obtain approproate binomial factor - accounts for samples
    binfac = facSort(nLin(i));
    
    if nLin(i) == 1 && any(fn.svec == tLin(i+1)) && singLinRem
        % Account for 1 lineage by just going to next time which would be a
        % sample time unless the simulation is at an end
        tsol = tLin(i+1);
        qsol = qev(i, :);
        disp(['Lineages prematurely fell to 1 at ' num2str(i)]);
    else
        % Solve linear ODEs continuously with setting of options, no Q
        if fn.id ~= 2
            %[tsol, qsol] = ode113(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, linBool),...
             %   [tLin(i) tLin(i+1)], qev(i, :)', options);
            [tsol, qsol] = ode23t(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, linBool),...
                [tLin(i) tLin(i+1)], qev(i, :)', options);
        else
            [tsol, qsol] = ode45(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, linBool),...
                [tLin(i) tLin(i+1)], qev(i, :)', options);
        end
    end
    
    % Normalise the posterior probabilities - only if linear ODE form
    if linBool
        for j = 1:size(qsol, 1)
            qsol(j, :) = qsol(j, :)/(sum(qsol(j, :)));
        end
    end
    
    % Assign the output values of time and posteriors
    qnoev(i, :) = qev(i, :);
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Calculate lamt and Nt across time explicitly
    lamtdiag{i} = zeros(length(tsol), m);
    for j = 1:length(tsol)
        % Calculate lamt across time explicitly
        [~, lamtdiag{i}(j, :)] = getTimeVaryingN(fn, tsol(j), binfac, xsetMx);
    end
    
    % Perturb the q posterior for the new event
    lampert = lamtdiag{i}(end, :);
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
    
    % Estimate the rate and population directly
    lamhatset{i} = zeros(1, length(tsol));
    for j = 1:length(tsol)
        % Parameters, rate and N(t)
        lamhatset{i}(j) = qsol(j, :)*lamtdiag{i}(j, :)';
    end
    Nhatset{i} = binfac./lamhatset{i};
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

% Calculate true values of population and rate across time
N = zeros(lenFull, 1);
lam = zeros(lenFull, 1);
for j = 1:lenFull
    % Calculate true N across time explicitly
    [N(j), ~] = getTimeVaryingN(fn, tn(j), binfac, x');
end
for i = 1:nData
    % Get lam across time
    if i < nData
        iddata = find(tn < tLin(i+1) & tn >= tLin(i));
    else
        % Accounts for last data point as no further coalescents
        iddata = find(tn <= tLin(i+1) & tn >= tLin(i));
    end
    lam(iddata) = facSort(nLin(i))./N(iddata);
end

% Calculate parameter estimates in general using xsetMx
xsetMx = fn.xsetMx;
% Conditional mean
xhat = qn*xsetMx';
xhatmean = xhat;
% Conditional variance
xhatvar = qn*(xsetMx'.^2) - xhat.^2;
xhatstd = sqrt(xhatvar);
xstdub = xhatmean + 2*xhatstd;
xstdlb = xhatmean - 2*xhatstd;
        
% Output structure
statsStruc.xhatmean = xhatmean;
statsStruc.xhatvar = xhatvar;
statsStruc.xhatstd = xhatstd;
statsStruc.xstdub = xstdub;
statsStruc.xstdlb = xstdlb;