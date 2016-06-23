% Function to simulate exp or con-exp-con functions with time rescaling
function [nLin, tLin, tcoal] = simExpPyb(svec, nvec, n, nData, x, isexp, facSort, nSampTimes)

% Simulate coalescent times for con-exp-con with rejection sampling and
% insertion of heterochronous sample times
tcoal = zeros(1, n); % defined so tcoal(1) = 0 = svec(1)
nLin = zeros(1, nData + 1);
tLin = zeros(1, nData + 1); % all lineage values
nLin(1) = nvec(1);
tLin(1) = svec(1);

% Count variables for simulation and booleans
iev = 1; % total no. loops
nLinCurr = nLin(1); % current no. lineages
tLinCurr = tLin(1); % start time
isamp = 1; % count no. samples <= nSampTimes
icoal = 0; % count no. coalescents <= n-1
sampVec = [svec inf]; % like svec but with inf for when exhaust samples in condition

% Define con-exp-con boundaries with modification for exp
if ~isexp
    texp1 = x(3);
    texp2 = x(3) + x(4);
    x5 = x(1)*exp(-x(2)*x(4));
else
    texp1 = 0;
    texp2 = inf;
    x5 = 0;
    x(3) = 0;
    x(4) = inf;
end

while(iev < nData+1)
    % Lineages fall to 1 before end then use next sample time and update
    if nLinCurr == 1 && isamp < nSampTimes + 1
        % Update to next sample
        sampTrue = 1;
    else
        % Get coalescent time for current no. lineages based on which
        % segment the function is in e.g. t <= x(3) => N(t) is x(1)
        if tLinCurr <= texp1 && ~isexp
            % First linear segment
            tLinCurr = tLinCurr + x(1)*exprnd(1)/facSort(nLinCurr);
        elseif tLinCurr >= texp2 && ~isexp
            % Second linear segment
            tLinCurr = tLinCurr + x5*exprnd(1)/facSort(nLinCurr);
        else
            % Exponential segment
            tLinCurr = x(3) + log(exp(x(2)*(tLinCurr - x(3))) + x(2)*x(1)*exprnd(1)/facSort(nLinCurr))/x(2);
        end
        % If next coalescent time is after next sample time then sample
        if tLinCurr > sampVec(isamp+1)
            % Update to next sample
            sampTrue = 1;
        else
            % A coalescent event has occurred
            nLinCurr = nLinCurr - 1;
            icoal = icoal + 1;
            sampTrue = 0;
            tcoal(icoal+1) = tLinCurr;
        end
    end
    
    % Sample update involves lineages and times from svec and nvec
    if sampTrue
        isamp = isamp + 1;
        nLinCurr = nLinCurr + nvec(isamp);
        tLinCurr = svec(isamp); % svec should never be exceeded
    end
    
    % Store data
    iev = iev + 1;
    nLin(iev) = nLinCurr;
    tLin(iev) = tLinCurr;
end

% Check simulation made use of all data
if(iev ~= isamp + icoal)
    error('Not all time data used in simulation');
end
% Check for consistency in event times
tcheck = sort([tcoal svec(2:end)]);
if ~all(tcheck == tLin)
    disp('The times are inconsistent');
end
if any(isnan(tLin))
    error('NaN values in the simulated times');
end