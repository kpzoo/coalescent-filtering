%% Function that focuses on con-exp-con functions sampled across time only
function [neffM, cap, lamM, lamV, est] = pybHeteroFn(n, nSampTimes, maxSamp, mi)

% Assumptions and Modifications
% - removed plots, Nhatz, xhatvar etc and marginalisations
% - difference from script is input n, nSampTimes and maxSamp
% - only pybus piecewise exponential demographic functions
% - introduces the effective lineage count idea
% - in heterochronous case divides samples uniformly %% <----- can change to other distributions later

% clc
% close all
% clearvars

% Set number of parameters and size of discretised grid
numRV = 4;
%mi = 15*ones(1, numRV);
if length(mi) ~= numRV
    error('mi is of incorrect length');
end
m = prod(mi);

% Set number of samples and assume uniformly divided in time
%n = 100;
%nSampTimes = 10;
nvec = round(n/nSampTimes)*ones(1, nSampTimes);
%maxSamp = 100;
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

% Set parameter space and true values but specify x and y-x vs y
minSpace = [10000 0.1 10 10];
maxSpace = [50000 0.5 50 50];

% Get the instance of the random variables for demographic function
xset = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    xset{i} = linspace(minSpace(i), maxSpace(i), mi(i));
    x(i) = datasample(xset{i}, 1);
end

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

% Define con-exp-con boundaries
texp1 = x(3);
texp2 = x(3) + x(4);
x5 = x(1)*exp(-x(2)*x(4));

while(iev < nData+1)
    % Lineages fall to 1 before end then use next sample time and update
    if nLinCurr == 1 && isamp < nSampTimes + 1
        % Update to next sample
        sampTrue = 1;
    else
        % Get coalescent time for current no. lineages based on which
        % segment the function is in e.g. t <= x(3) => N(t) is x(1)
        if tLinCurr <= texp1
            % First linear segment
            tLinCurr = tLinCurr + x(1)*exprnd(1)/facSort(nLinCurr);
        elseif tLinCurr >= texp2
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

% Get N(t) con-exp-con in time - true function 
tz = linspace(tcoal(1), tcoal(end), 1000);
Nz = zeros(size(tz));
for i = 1:length(tz)
    I1 = tz(i) <= x(3);
    I2 = tz(i) > x(3) & tz(i) < x(4) + x(3);
    I3 = tz(i) >= x(4) + x(3);
    Nz(i) = x(1).*I1 + x(1).*exp(-x(2).*(tz(i) - x(3))).*I2 + ...
        x(1).*exp(-x(2).*x(4)).*I3;
end

% Check for case when x5 goes to zero and tcoal is not distinct
if any(diff(tcoal) == 0)
    error('The coalescent times are not distinct');
end

% Effective no. lineages in time with tz as time and nz the lineage number
nz = zeros(size(tz));
nz(1) = nLin(1);
for i = 2:length(tz)
    nz(i) = nLin(find(tLin <= tz(i), 1, 'last'));
end

% Effective number of lineages and coalescent rate and capacity - these
% indices may become more useful in batch testing different schemes
neffz = 0.5 + sqrt(0.25 + nz.*(nz - 1)./Nz);
neffM = trapz(tz, neffz)/range(tz);
lamz = nz.*(nz - 1)./Nz/2;
lamM = trapz(tz, lamz)/range(tz);
lamV = trapz(tz, lamz.*lamz)/range(tz) - lamM^2;
cap = lamM*log(max(lamz)/lamM);

%% Snyder Filtering of the coalescent-heterochronous data

% Initialise with uniform prior
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
%lamtdiag = cell(1, 1);

% Time saving calculations
cond1 = xsetMx(3, :);
cond2 = xsetMx(4, :) + xsetMx(3, :);

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
         [tsol, qsol] = ode113(@(ts, y) odeSnyPyb(ts, y, binfac, xsetMx),...
             [tLin(i) tLin(i+1)], qev(i, :)', options);
        
        % Assign the output values for time and posteriors
        qnoev(i, :) = qev(i, :);
        qset{i} = qsol;
        tset{i} = tsol;
        elemLen(i) = length(tsol);
        
%         % Calculate lamt and Nt across time explicitly
%         lamtdiag{i} = zeros(length(tsol), m);
%         for j = 1:length(tsol)
%             % Calculate lamt across time explicitly
%             I1 = tsol(j) <= cond1;
%             I2 = tsol(j) > cond1 & tsol(j) < cond2;
%             I3 = tsol(j) >= cond2;
%             Ntemp = xsetMx(1, :).*I1 + xsetMx(1, :).*exp(-xsetMx(2, :).*(tsol(j) - xsetMx(3, :))).*I2 + ...
%                 xsetMx(1, :).*exp(-xsetMx(2, :).*xsetMx(4, :)).*I3;
%             lamtdiag{i}(j, :) = binfac./Ntemp;
%         end
        
        % Only perturb if there is a new event as opposed to sample
        if ~any(tLin(i+1) == svec)
            % Perturb the q posterior for the new event %%% <-----------------maybe lampert needs (i+1)th binfac
            %lampert = lamtdiag{i}(end, :);
            
            % Only calculate rate diagonal when need to perturb
            I1 = tsol(end) <= cond1;
            I2 = tsol(end) > cond1 & tsol(end) < cond2;
            I3 = tsol(end) >= cond2;
            Ntemp = xsetMx(1, :).*I1 + xsetMx(1, :).*exp(-xsetMx(2, :).*(tsol(end) - xsetMx(3, :))).*I2 + ...
                xsetMx(1, :).*exp(-xsetMx(2, :).*xsetMx(4, :)).*I3;
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

% Conditional mean
xhat = qn*xsetMx';
% xhatmean = xhat;
% % Conditional variance
% xhatvar = max(0, qn*(xsetMx'.^2) - xhat.^2); %% negative values have no meaning
% xhatstd = sqrt(xhatvar);
% xstdub = xhatmean + 2*xhatstd;
% xstdlb = xhatmean - 2*xhatstd;

% Display structure for parameters
est.x = x;
est.xhat = xhat(end, :);



%% Post calculations

% % Marginalise the posterior across sim times - take 20 of them
% nPoster = 20;
% idpost = round(linspace(1, length(tn), nPoster));
% qmarg = cell(1, nPoster);
% probSums = cell(1, nPoster);
% for i = 1:nPoster
%     [qmarg{i}, probSums{i}] = marginalise(numRV, IDMx, qn(idpost(i), :), mi);
% end
% qnlast = qn(end, :);
% 
% % Reorganise into probability matrices for each variable
% qnp = cell(1, numRV);
% for i = 1:numRV
%     qnp{i} = zeros(nPoster, mi(i));
%     for j = 1:nPoster
%         qnp{i}(j, :) = qmarg{j}{i};
%     end
% end
% tnmin = min(tn);
% tnmax = max(tn);
% tnlen = length(tn);
% 
% % Final estimates - only works for inhomogeneous MCs
% fn.xsetMx = xsetMx;
% fn.id = 6;
% [Nzhat, Nbnd] = getPopBackward(tz, qnlast, fn);
% 
% % Plots to visualise results
% 
% % Lineage through time plot with sample times
% figure;
% stairs(tLin, nLin);
% hold on
% plot(repmat(svec, [2 1]), [zeros(size(svec)); n*ones(size(svec))], 'r--', 'linewidth', 2);
% hold off
% legend('lineages', 'sample times', 'location', 'best');
% xlabel('time');
% ylabel('no. lineages');
% title(['Lineages through time with no. sample times = ' num2str(nSampTimes)]);
% xlim([tLin(1) tLin(end)]);
% ylim([min(nLin) max(nLin)]);
% 
% 
% % Posteriors and final estimates
% figure;
% for i = 1:numRV
%     subplot(numRV, 1, i);
%     plot(xset{i}, qnp{i}(1, :), 'b-', xset{i}, qnp{i}(end, :), 'ko-',...
%         xset{i}, qnp{i}(2:end-1, :), 'g--');
%     hold on
%     plot([xhat(end, i) xhat(end, i)], [0 max(max(qnp{i}))], 'k');
%     hold off
%     % Condition to account for variables with only one value
%     if min(xset{i}) ~= max(xset{i})
%         xlim([min(xset{i}) max(xset{i})]);
%     end
%     xlabel(['x_' num2str(i)]);
%     ylabel(['P(x_' num2str(i) '|data)']);
%     legend('prior', 'posterior', 'intermediates', 'location', 'best');
%     title(['Evolution of posteriors: param est = ' num2str(xhat(end, i))]);
% end
% 
% 
% % The parameters estimated and the std deviation bounds
% figure;
% for i = 1:numRV
%     subplot(numRV, 1, i);
%     plot(tn, x(i)*ones(size(tn)), 'r', tn, xhatmean(:, i), 'k', tn, xstdlb(:, i),...
%         'b', tn, xstdub(:, i), 'b');
%     xlabel('time');
%     ylabel(['x_' num2str(i) ' in N(t) =' 'x_1e^-x_2t']);
%     title(['Estimate of x_' num2str(i) 'for [n m] = [' num2str(n) ' ' num2str(m) ']']);
%     legend('true', 'estimate', '-2 std', '+2 std', 'location', 'best');
%     xlim([tnmin tnmax]);
% end
% 
% % Final N(t) estimate from final parameters with coalescent points
% vx_coal = repmat(tcoal, 2, 1);
% v_coal = [zeros(1, n); max(max(Nbnd))*ones(1, n)];
% vx_samp = repmat(svec(2:end), 2, 1);
% v_samp = [zeros(1, nSampTimes-1); max(max(Nbnd))*ones(1, nSampTimes-1)];
% figure;
% plot(tz, Nz, 'r', 'LineWidth', 2);
% hold on
% plot(tz, Nzhat, 'k', 'LineWidth', 2);
% plot(tz, Nbnd, 'g--', 'LineWidth', 2);
% plot(vx_samp, v_samp, 'm');
% plot(vx_coal, v_coal, 'k:');
% hold off
% legend('N(t)', 'Nhat(t)', 'lower bound', 'upper bound', 'sample times', 'coalescents', 'location', 'best');
% xlim([tcoal(1) tcoal(end)]);
% xlabel('time');
% ylabel('N(t) estimate and bounds');
% title('Population estimate using final posterior');
% grid;
