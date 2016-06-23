% Perform Snyder filtering for constant N

% Assumptions
% - the effective population is constant and can take m random values
% - the data is k-coalescent times and the process self-exciting

clc 
close all
clear all

%% Filter the simulated event times from the appropriate Poisson process

% Set initial parameters: q0 is prior, N is true value, m is dimension of N
% space and n is no. samples
m = 10;
q0 = ones(1, m)/m;
n = 1000;
nData = n-1;

% Get the instance of random variable N and initial rate set vector without
% the self excitation, lvec
%Nset = 50*(1:m) + 100;
Nset = n*[100:200:2000];
N = randsample(Nset, 1);
lvec = 1./Nset;

% Get rates of coalescent and the waiting times, twait as well as
% cumulative times and the MRCA time
nset = n:-1:2;
rates = (nset.*(nset-1))/(2*N);
twait = exprnd(1./rates);
tcoal = cumsum([0 twait]);
tmrca = tcoal(end);

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nData+1, m);
qnoev = qev;
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
Nhatset = cell(1, 1);
options = odeset('NonNegative', 1:m);
elemLen = zeros(1, nData);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain rate matrix, lam
    fac = 0.5*(n - i + 1)*(n - i);
    lam = diag(fac*lvec);
   
    % Solve linear ODEs continuously with setting of options, no Q
    [tsol, qsol] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam),...
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
    
    % Perturb the q posterior for the new event
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :)*lam./sum(qev(i+1, :)*lam);
    
    % Estimate the population parameter N from the posterior
    Nhatset{i} = fac./(qsol*diag(lam));
end


% Get full length of ODE solution data and assign appending vectors
lenFull = sum(elemLen);
stop = 0;
qn = -ones(lenFull, m);
tn = -ones(lenFull, 1);
Nhat = -ones(lenFull, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    tn(start:stop) = tset{i};
    qn(start:stop, :) = qset{i};
    Nhat(start:stop, :) = Nhatset{i};
end

%% Convert data in MSE estimates

% Resample the time series by first removing any initial times and deciding
% a sample interval, delt and get raw error, en
tn = tn - tn(1);
delt = mean(diff(tn))/20;
nSamps = range(tn)/delt;
en = N - Nhat;

% Remove duplicate points - the +1 in the id chooses the latter of point
dtn = diff(tn);
id = find(dtn == 0) + 1; 
lentn = length(tn);

% Obtain altered data vectors with duplicate time indices removed
idset = 1:lentn;
trunc = setdiff(idset', id);
tm = tn(trunc);
em = en(trunc);

% Separate samples in sets of Tsetsize
Tsetsize = 20000;
nSetSamps = floor(nSamps/Tsetsize);
tSetSamps = range(tn)/nSetSamps;

% Linear interpolation variable declaration for en, tlim(1) set for no bias
sum_em = 0;
sum_emSq = 0;
sum_len = 0;
tlim = -ones(1, nSetSamps+1);
tlim(1) = 0;

% Loop across sets iteratively calculating statistics
for i = 2:nSetSamps+1
    % Define time limits and obtain relevant section of data
    tlim(i) = tSetSamps*(i-1);
    idtemp = find(tm >= tlim(i-1) & tm < tlim(i));
    ttemp = tm(idtemp);
    etemp = em(idtemp);
    
    % Obtain sample times and interpolate em data to the sampled times
    tsamp = ttemp(1):delt:ttemp(end);
    esamp = nakeinterp1(ttemp, etemp, tsamp');
    
    lenSamp = length(esamp);
    if lenSamp ~= length(tsamp)
        assignin('base', 'tsamp', tsamp);
        assignin('base', 'eSamp', esamp);
        error(['The sampling of the error curve failed at i = ' num2str(i)]);
    end
    
    % Obtain iterative sums for statistics
    esampSq = esamp.*esamp;
    sum_em = sum_em + sum(esamp);
    sum_emSq = sum_emSq + sum(esampSq);
    sum_len = sum_len + lenSamp;
end

% Check that the correct number of points were obtained and calculate means
if sum_len > nSamps 
    assignin('base', 'tlim', tlim);
    error(['Inconsistent sample size: [nSamps sum_len] = ' [num2str(nSamps)...
        ' ' num2str(sum_len)]]);
else
    em_mean = sum_em/sum_len;
    em_mse = sum_emSq/sum_len;
    em_var = em_mse - em_mean^2;
end

% Assign stats set with indicator string to indicate data order
statSet.order = {'mean', 'var', 'mse'};
statSet.val = [em_mean em_var em_mse];

%% Plotting results

% Plot the estimate of N
figure;
plot(tn, Nhat, tn, N*ones(size(tn)), 'k');
ylim([min(min(Nhat), N) - 0.1*N max(max(Nhat), N) + 0.1*N]);
legend('Nhat', 'N', 'location', 'best')
xlabel('time');
ylabel('effective population');

% Plot the absolute percentage error
Sm = 100*abs(em)/N;
figure;
plot(tm, Sm, 'o');
xlabel('time');
ylabel('perc |N - Nhat|');