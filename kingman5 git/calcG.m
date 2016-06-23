% Calculate a series of test parameters for coalescents
clear all
clc
close all

%% Calculate G from the capacity equation
n = 2:100;
lenn = length(n);
G = zeros(1, lenn);

for i = 1:lenn
    jmax = n(i);
    lenj = jmax-1;
    ff = zeros(1, lenj);
    for j = 1:lenj
        ff(j) = nchoosek(jmax - (j-1), 2);
    end
    G(i) = sum(1./ff);
end

%% Plot capacity as a function of n, assuming G = 2
n = 2:1000;
N = 10^6;
C = (n-1).*log(n/2)/(2*N);
dC = (1 + log(n/2) - 1./n)/(2*N);


%% Fano inequality for error

% Define sample no. n and message no. m
n = 10:100;
m = 2;
theta = 10; % if m > 2 then need to define theta1, theta2 etc

% Euler constant and harmonic sum
g = 0.577215;
harsum = size(n);
for i = 1:length(n)
    div = 1:n(i);
    div = 1./div;
    harsum(i) = sum(div);
end

% Only for m = 2 case

%Hxy = (n-1).*(g + 1./n - harsum) + n + log(factorial(n-1)/theta);
Hxy = (log(gamma(n)/theta) + n - (n-1).*psi(n))/m;
%Perr = (Hxy - 1)/log(m);
p1 = 0:0.01:1;
Hp1 = p1.*log(1./p1);

figure;
plot(Hxy);
hold
plot(Hp1);
hold

