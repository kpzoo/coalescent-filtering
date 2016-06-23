% Function to perform Snyder filtering for con-exp-con and exponential N(t)
function [xhat, xhatmean, xstdlb, xstdub, qnlast] = pybexpFilterDel(mi, nData, isexp, xsetMx, tLin, nLin, facSort, svec, delRV)

% Assumptions and modifications
% - added allowance for fixed x(3) and x(4)
% - maybe lampert needs (i+1)th binfac <---------------------

% Initialise with uniform prior of appropriate dimension
if isexp
    % Exponential is bivariate
    m = prod(mi(1:2));
    xsetMx = xsetMx(1:2, 1:m);
else
    % Con-exp-con has 4 parameters
    m = prod(mi);
end
q0 = ones(1, m)/m;

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nData + 1, m);
qnoev = qev;
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
%options = odeset('NonNegative', 1:m);
options= odeset('RelTol', 1e-9, 'NonNegative', 1:m);
elemLen = zeros(1, nData);

% Time saving calculations
if ~isexp
    if delRV.bool
        cond1 = delRV.x3;
        cond2 = delRV.x3 + delRV.x4;
    else
        cond1 = xsetMx(3, :);
        cond2 = xsetMx(4, :) + xsetMx(3, :);
    end
else
    cond1 = -1;
    cond2 = -1;
end

% Run Snyder filter across the time series accounting for 
for i = 1:nData
    
    % Obtain appropriate binomial factor - accounts for samples
    binfac = facSort(nLin(i));
    
    % Solve posterior if nLin >= 2 else do nothing
    if binfac == 0
        % Maintain posterior at last value if lineages fall to 1
        qev(i+1, :) = qev(i, :);
        disp(['Lineages prematurely fell to 1 at ' num2str(i)]);
    else
         % Solve Snyder ODEs for posterior
         if ~isexp
             % ODE integration con-exp-con
             if delRV.bool
                 [tsol, qsol] = ode113(@(ts, y) odeSnyPybDel(ts, y, binfac, xsetMx, delRV),...
                     [tLin(i) tLin(i+1)], qev(i, :)', options);
             else
                 [tsol, qsol] = ode113(@(ts, y) odeSnyPyb(ts, y, binfac, xsetMx),...
                     [tLin(i) tLin(i+1)], qev(i, :)', options);
             end
         else
             % ODE integration exponential
             [tsol, qsol] = ode113(@(ts, y) odeSnyExp(ts, y, binfac, xsetMx),...
                 [tLin(i) tLin(i+1)], qev(i, :)', options);
         end
        
        % Assign the output values for time and posteriors
        qnoev(i, :) = qev(i, :);
        qset{i} = qsol;
        tset{i} = tsol;
        elemLen(i) = length(tsol);
        
        % Only perturb if there is a new event as opposed to sample
        if ~any(tLin(i+1) == svec)
            
            % Calculate rate diagonal when need to perturb q for new event
            if ~isexp
                % Con-exp-con function
                I1 = tsol(end) <= cond1;
                I2 = tsol(end) > cond1 & tsol(end) < cond2;
                I3 = tsol(end) >= cond2;
                if ~delRV.bool
                    Ntemp = xsetMx(1, :).*I1 + xsetMx(1, :).*exp(-xsetMx(2, :).*(tsol(end) - xsetMx(3, :))).*I2 + ...
                        xsetMx(1, :).*exp(-xsetMx(2, :).*xsetMx(4, :)).*I3;
                else
                    % Removed 2 variables
                    Ntemp = xsetMx(1, :).*I1 + xsetMx(1, :).*exp(-xsetMx(2, :).*(tsol(end) - delRV.x3)).*I2 + ...
                        xsetMx(1, :).*exp(-xsetMx(2, :).*delRV.x4).*I3;
                end
            else
                % Exponential function
                Ntemp = (xsetMx(1, :).*exp(-xsetMx(2, :)*tsol(end)));
            end
            lampert = binfac./Ntemp;
            
            qev(i+1, :) = qsol(end, :);
            tev(i+1) = tsol(end);
            qev(i+1, :) = qev(i+1, :).*lampert./(qev(i+1, :)*lampert');
        else
            % Maintain posterior since new samples only
           qev(i+1, :) = qsol(end, :);
        end
        
    end
end

% Get full length of ODE solution data and assign appending vectors
lenFull = sum(elemLen);
stop = 0;
qn = -ones(lenFull, m);
tn = -ones(lenFull, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    tn(start:stop) = tset{i};
    qn(start:stop, :) = qset{i};
end

% Conditional mean and last posterior
xhat = qn*xsetMx';
xhatmean = xhat;
qnlast = qn(end, :);
% Conditional variance
xhatvar = max(0, qn*(xsetMx'.^2) - xhat.^2); %% negative values have no meaning
xhatstd = sqrt(xhatvar);
xstdub = xhatmean + 2*xhatstd;
xstdlb = xhatmean - 2*xhatstd;