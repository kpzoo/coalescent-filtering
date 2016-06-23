% Function to perform Snyder filtering for heterochronous data
function [xhat, Nhat, lamhat, lam, N, tn, statsStruc, qn] = snyderHetero2(fn, q0, events)

% Assumptions and modifications
% - nData is no longer n-1 but includes sample time events
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

% Extract heterochronous inputs
tLin = events.tLin;
nLin = events.nLin;

% Data here is coalescent and sample events from second sample time on
nD = length(tLin)-1;
tsamp = fn.svec(2:end);

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nD+1, m);
qnoev = qev;
qev(1, :) = q0;
tev = zeros(nD, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
lamhatset = cell(1, 1);
Nhatset = cell(1, 1);
options = odeset('NonNegative', 1:m);
%options= odeset('RelTol', 1e-9, 'NonNegative', 1:m);
elemLen = zeros(1, nD);
lamtdiag = cell(1, 1);

% Run Snyder filter across the time series accounting for 
for i = 1:nD
    
    % Obtain appropriate binomial factor - accounts for samples
    binfac = facSort(nLin(i));
    
    % Solve posterior if nLin >= 2 else do nothing
    if binfac == 0
        % Maintain posterior at last value if lineages fall to 1
        qev(i+1, :) = qev(i, :);
        disp(['Lineages prematurely fell to 1 at ' num2str(i)]);
    else
         % Solve Snyder ODEs for posterior
        if fn.id ~= 2
            [tsol, qsol] = ode113(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, 0),...
                [tLin(i) tLin(i+1)], qev(i, :)', options);
        else
            [tsol, qsol] = ode45(@(ts, y) odeSnyder2(ts, y, binfac, fn, xsetMx, 0),...
                [tLin(i) tLin(i+1)], qev(i, :)', options);
        end
        
        % Assign the output values for time and posteriors
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
        
        % Only perturb if there is a new event as opposed to sample
        if ~any(tLin(i+1) == tsamp)
            % Perturb the q posterior for the new event %%% <-----------------maybe lampert needs (i+1)th binfac
            lampert = lamtdiag{i}(end, :);
            %[~, lampert] = getTimeVaryingN(fn, tsol(end), facSort(nLin(i)), xsetMx);
            qev(i+1, :) = qsol(end, :);
            tev(i+1) = tsol(end);
            qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
        else
            % Maintain posterior since new samples only
           qev(i+1, :) = qsol(end, :);
        end
        
        % Estimate the rate and population directly
        lamhatset{i} = zeros(1, length(tsol));
        for j = 1:length(tsol)
            % Parameters, rate and N(t)
            lamhatset{i}(j) = qsol(j, :)*lamtdiag{i}(j, :)';
        end
        Nhatset{i} = binfac./lamhatset{i};
        
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
for i = 1:nD
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
for i = 1:nD
    % Get lam across time
    if i < nD
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
xhatvar = max(0, qn*(xsetMx'.^2) - xhat.^2); %% negative values have no meaning
xhatstd = sqrt(xhatvar);
xstdub = xhatmean + 2*xhatstd;
xstdlb = xhatmean - 2*xhatstd;
        
% Output structure
statsStruc.xhatmean = xhatmean;
statsStruc.xhatvar = xhatvar;
statsStruc.xhatstd = xhatstd;
statsStruc.xstdub = xstdub;
statsStruc.xstdlb = xstdlb;