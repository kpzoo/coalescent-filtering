% Perform Snyder filtering for constant N

% Assumptions and Modifications
% - the effective population is constant and can take m random values
% - the data is k-coalescent times and the process self-exciting
% - added option to scale data and use a constant lam and compare

clc 
close all
clear all

%% Filter the simulated event times from the appropriate Poisson process

% Set initial parameters: q0 is prior, N is true value, m is dimension of N
% space and n is no. samples
m = 100;
q0 = ones(1, m)/m;
n = 100;
nData = n-1;

% Get the instance of random variable N and initial rate set vector without
% the self excitation, lvec
%Nset = 50*(1:m) + 100;
Nset = n*linspace(100, 1000, m);
N = randsample(Nset, 1);
lvec = 1./Nset;

% Get rates of coalescent and the waiting times, twait as well as
% cumulative times and the MRCA time
nset = n:-1:2;
fac = nset.*(nset-1)/2;
rates = fac/N;
twait = exprnd(1./rates);
tcoal = cumsum([0 twait]);
tmrca = tcoal(end);

% Scale the data for comparative filtering and obtain lam for it
twait2 = twait.*fac;
%twait2 = exprnd(N, size(twait)); % <--------- Test with t2 data
tcoal2 = cumsum([0 twait2]);
lam2 = diag(lvec);

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nData+1, m);
qnoev = qev;
qev(1, :) = q0;
tev = zeros(nData, 1);

% Duplicate vectors for scaled data
qev2 = zeros(nData+1, m);
qnoev2 = qev2;
qev2(1, :) = q0;
tev2 = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);
Nhatset = cell(1, 1);
qset2 = cell(1, 1);
tset2 = cell(1, 1);
Nhatset2 = cell(1, 1);
options = odeset('NonNegative', 1:m);
elemLen = zeros(1, nData);
elemLen2 = zeros(1, nData);
Nev = zeros(1, nData);
Nev2 = zeros(1, nData);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain rate matrix, lam with account for scaling
    lam = diag(fac(i)*lvec);
   
    % Solve linear ODEs continuously with setting of options, no Q
    % Solution with scaled data
    [tsol2, qsol2] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam2),...
        [tcoal2(i) tcoal2(i+1)], qev2(i, :)', options);
    % Solution with direct data
    [tsol, qsol] = ode113(@(ts, y) odeSnyLin(ts, y, [], lam),...
        [tcoal(i) tcoal(i+1)], qev(i, :)', options);
    
    % Normalise the posterior probabilities
    for j = 1:size(qsol, 1)
        qsol(j, :) = qsol(j, :)/(sum(qsol(j, :)));
    end
    for j = 1:size(qsol2, 1)
        qsol2(j, :) = qsol2(j, :)/(sum(qsol2(j, :)));
    end
    
    % Assign the output values of time and posteriors
    qnoev(i, :) = qev(i, :);
    qset{i} = qsol;
    tset{i} = tsol;
    elemLen(i) = length(tsol);
    
    % Assign the output values of time and posteriors for scaled data
    qnoev2(i, :) = qev2(i, :);
    qset2{i} = qsol2;
    tset2{i} = tsol2;
    elemLen2(i) = length(tsol2);
    
    % Perturb the q posterior for the new event
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :)*lam./sum(qev(i+1, :)*lam);

    % Scaled data updates
    qev2(i+1, :) = qsol2(end, :);
    tev2(i+1) = tsol2(end);
    qev2(i+1, :) = qev2(i+1, :)*lam2./sum(qev2(i+1, :)*lam2);
    
    % Estimate the population parameter N from the posterior
    Nhatset{i} = fac(i)./(qsol*diag(lam));
    Nhatset2{i} = 1./(qsol2*diag(lam2));
    Nev(i) = fac(i)./(qev(i+1, :)*diag(lam));
    Nev2(i) = 1./(qev2(i+1, :)*diag(lam2));
end


% Get full length of ODE solution data and assign appending vectors
lenFull = sum(elemLen);
lenFull2 = sum(elemLen2);
stop = 0;
stop2 = 0;
qn = -ones(lenFull, m);
tn = -ones(lenFull, 1);
Nhat = -ones(lenFull, 1);
qn2 = -ones(lenFull2, m);
tn2 = -ones(lenFull2, 1);
Nhat2 = -ones(lenFull2, 1);

% Append the cell based posterior and time data into a single structure
for i = 1:nData
    % Loop calculates the start and end points along the array successively
    % and then assigns the appropriate cell element
    start = stop + 1;
    stop = elemLen(i) + start - 1;
    tn(start:stop) = tset{i};
    qn(start:stop, :) = qset{i};
    Nhat(start:stop, :) = Nhatset{i};
    
    % Append data sets for scaled set
    start2 = stop2 + 1;
    stop2 = elemLen2(i) + start2 - 1;
    tn2(start2:stop2) = tset2{i};
    qn2(start2:stop2, :) = qset2{i};
    Nhat2(start2:stop2, :) = Nhatset2{i};
end

%% Convert data in MSE estimates and plot comparisons

% MSE calculation and time scales
[statSet, em, tm] = coalMSE(tn, N, Nhat);
[statSet2, em2, tm2] = coalMSE(tn2, N, Nhat2);
tnmin = min([tn; tn2]);
tnmax = max([tn; tn2]);

% Plot the estimates of N
figure;
plot(tn, Nhat, tn2, Nhat2, [tnmin; tnmax], [N; N], 'k');
ylim([min(min([Nhat; Nhat2]), N) - 0.1*N max(max([Nhat; Nhat2]), N) + 0.1*N]);
legend('Nhat', 'Nhat2', 'N', 'location', 'best')
xlabel('time');
ylabel('effective population');

% Plot the absolute percentage error
Sm = 100*abs(em)/N;
Sm2 = 100*abs(em2)/N;
figure;
plot(tm, Sm, 'bo', tm2, Sm2, 'ro');
xlabel('time');
ylabel('perc |N - Nhat|');
legend('unscaled', 'scaled', 'location', 'best')

% Check against classical MMSE estimators = suff stat/(ndata - 2)
T1 = sum(fac.*twait);
T2 = sum(twait2);
Nclass = T1/(nData - 2);
Nclass2 = T2/(nData - 2);
Nsny = Nhat(end);
Nsny2 = Nhat2(end);

% Display comparisons
disp(['N is ' num2str(N)]);
disp(['Nhat(end) is ' num2str(Nsny)]);
disp(['Nhat2(end) is ' num2str(Nsny2)]);
disp(['Nclass is ' num2str(Nclass)]);
disp(['Nclass2 is ' num2str(Nclass2)]);

% Get squared error ratios
RsqErr = [(1 - Nsny/N)^2 (1 - Nsny2/N)^2 (1 - Nclass/N)^2 (1 - Nclass2/N)^2];
disp(['Squared error ratios are: ' num2str(RsqErr)]);

% Get cumulative classical estimators
t1 = cumsum(fac.*twait);
t2 = cumsum(twait2);
rr = (1:nData) - 2;
Nc = t1./rr;
Nc2 = t2./rr;

% Remove the division by zero and negative values and make Nevs consistent
Nc = Nc(3:end);
Nc2 = Nc2(3:end);
Nevf = Nev(3:end);
Nev2f = Nev2(3:end);
newlen = length(Nc);

% Get MSE ratio of cumulative classical and Snyder estimators
mseCum = zeros(4, newlen);
mseCum(1, :) = (1 - Nevf/N).^2;
mseCum(2, :) = (1 - Nev2f/N).^2;
mseCum(3, :) = (1 - Nc/N).^2;
mseCum(4, :) = (1 - Nc2/N).^2;

%%% NOTE - Nhat and Nhat2 appear to have the same values but across
%%% different time scales
figure;
plot(Nhat);
hold on
plot(Nhat2, 'r');
hold off

% Plot estimators across events or data samples
figure;
semilogy(1:newlen, mseCum(1, :), 'b', 1:newlen, mseCum(3, :), 'r');
xlabel('events');
ylabel('MSE estimates');
legend('Nev', 'Nc', 'location', 'best');