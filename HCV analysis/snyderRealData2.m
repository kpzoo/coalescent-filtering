% Function to perform Snyder filtering and estimation
function [xhat, Nhat, lamhat, tn, nLin, statsStruc, qn] = snyderRealData2(fn, q0, tcoal, linBool)

% Assumptions and modifications
% - added non-linear snyder form
% - added the tcoal(i) == tcoal(i+1) case
% - modified version of snyder to account for real data
% - removed outputs of N and lam as not known as well as 2 RV special calc
% - works for multivariable functions defined in fn
% - includes new functions for calculating lamdiag and N
% - uses matrices xsetMx etc defined in fn

% Extract inputs
fac = fn.fac;
nData = fn.nData;
m = fn.m;
xsetMx = fn.xsetMx;

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

% Cells for true values at evaluated times
lamtdiag = cell(1, 1);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain approproate binomial factor
    binfac = fac(i);
    
    % Solve linear ODEs continuously with setting of options, no Q
    if tcoal(i) ~= tcoal(i+1)
        if fn.id ~= 2
            [tsol, qsol] = ode113(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, linBool),...
                [tcoal(i) tcoal(i+1)], qev(i, :)', options);
        else
            [tsol, qsol] = ode45(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, linBool),...
                [tcoal(i) tcoal(i+1)], qev(i, :)', options);
        end
    else
        % Identical coalescent times - just assign last value
        tsol = tcoal(i+1);
        qsol = qev(i, :);
        disp(['Duplicate coalescent at ' num2str(i)]);
    end
    
    % Normalise the posterior probabilities
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
    
    % Calculate lamt across time explicitly;
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

% Calculate number of lineages through time and binomial factors 
nLin = zeros(lenFull, 1);
facLin = zeros(lenFull, 1);
for i = 1:nData
    % Get binfac and no. lineages across time and so lam from N
    if i < nData
        iddata = find(tn < tcoal(i+1) & tn >= tcoal(i));
    else
        % Accounts for last data point as no further coalescents
        iddata = find(tn <= tcoal(i+1) & tn >= tcoal(i));
    end
    facLin(iddata) = fac(i)*ones(size(iddata));
    nLin(iddata) = (nData+1 - i)*ones(size(iddata));
end

% Calculate parameter estimates - conditional mean and variance
xhat = qn*xsetMx';
xhatmean = xhat;
xhatvar = qn*(xsetMx'.^2) - xhat.^2;
xhatstd = sqrt(xhatvar);
xstdub = xhatmean + 2*xhatstd;
xstdlb = xhatmean - 2*xhatstd;
    
% Output structure
statsStruc.facLin = facLin;
statsStruc.xhatmean = xhatmean;
statsStruc.xhatvar = xhatvar;
statsStruc.xhatstd = xhatstd;
statsStruc.xstdub = xstdub;
statsStruc.xstdlb = xstdlb;