% Test code to determine MSE convergence of two types of estimators
clear all
clc
close all

% Note that n = k samples means coalescent times from k+1 to 2 lineages

% Set a number of samples
n = 3:5000;
lenn = length(n);

% Simulate the conditioning data for the maximum length and keep N >> nmax
nmax = n(end);
N = 100*nmax;
theta = 1/N;

% First data set is iid from first coalescent times
D1 = exprnd(1/theta, [1, nmax]);
% Second data set is from up to n coalescents
D2 = zeros(1, nmax);
fac = D2;
for i = 1:nmax
    fac(i) = nchoosek(i+1, 2);
    D2(i) = exprnd(1/(fac(i)*theta));
end

% Calculate the sufficient statistics
T1 = zeros(1, lenn);
T2 = T1;

% Obtain sufficient statistics for each sample length
for i = 1:lenn
    T1(i) = sum(D1(1:n(i)));
    T2(i) = sum(fac(1:n(i)).*D2(1:n(i)));
end

% Get the estimates with sample length <------------- put n-2 for MMSE
a = 2;
N1 = T1./(n-a);
N2 = T2./(n-a);

% Calculate abs percentage error and MSE
perc1 = 100*abs(N*ones(size(n)) - N1)/N;
perc2 = 100*abs(N*ones(size(n)) - N2)/N;
mse1 = (N*ones(size(n)) - N1).^2;
mse2 = (N*ones(size(n)) - N2).^2;

% Plot estimates
figure;
plot(n, N1, n, N2, n, N*ones(size(n)), 'k');
legend('location', 'best', 't_2 data', 't_n data', 'true value');
xlabel('no. samples');
ylabel('N and Nhat');

figure;
semilogy(n, mse1, n, mse2);
xlabel('no. samples');
ylabel('MSE');
legend('location', 'best', 't_2 data', 't_n data');

figure;
semilogy(n, perc1, n, perc2);
xlabel('no. samples');
ylabel('percentage error');
legend('location', 'best', 't_2 data', 't_n data');

%%%%%%%%%% Check %%%%%%%%%%
% Does the coalescent time data involve waiting times or sums of waiting
% times? Ans the times between else the likelihood analysis collapses as
% the independence criterion is lost - but Snyder will only see sum of
% times really - thought the estimator updates on events so maybe not