% Script to pull out many garli coalescents and then get estimates in batch
clc 
close all
clearvars

% Assumptions and modifications
% - assumes fnid = 6 as piecewise exponential on HCV data

% Define file name (text file) and extract the trees in newick form
name = '68set.text';
[nTree, treeStr] = getNexusStrings(name);

% Control booleans and variables
bool.plotIndiv = 0;
bool.remDuplicates = 1;
bool.mu = 0.79*(10^-3);
bool.fnid = 6;

% Number of time values
ntz = 1000;
t = zeros(nTree, ntz);
Nmean = zeros(nTree, ntz);

% Variables to store coalescents, estimates and populations
nold = zeros(1, nTree);
n = zeros(1, nTree);
tIter = zeros(1, nTree);
tcoal = cell(1, nTree);
xest = n;
yest = n;
NAest = n;
NCest = n;
rest = n;
TMRCAest = n;

% Loop across trees, obtain coalescent times and perform inference
for ii = 1:nTree
    % Main code for filtering and estimation
    nameNewick = treeStr{ii};
    [tcoal{ii}, est, tIter(ii), Nmean(ii, :), t(ii, :), nold(ii), n(ii)] = multiFilterGarli(nameNewick, bool, ntz);
    % Parameter estimates
    xest(ii) = est.x(2);
    yest(ii) = est.y(2);
    NAest(ii) = est.NA;
    NCest(ii) = est.NC(2);
    rest(ii) = est.r(2);
    TMRCAest(ii) = est.TMRCA;
end

% Parameter estimates from Pybus2003
if all(nold == 68)
    % Data from 68 sequences
    oli.NC = [4095 10310 18960];
    oli.r = [0.075 0.264 0.620];
    oli.x = [1941 1953 1966];
    oli.y = [1924 1934 1943];
    oli.NA = [153 245 345];
    oli.TMRCA = [1258 1374 1481];
elseif all(nold == 63)
    % Data from 63 sequences
    oli.NC = [3323 8779 15780];
    oli.r = [0.072 0.237 0.564];
    oli.x = [1941 1953 1966];
    oli.y = [1922 1932 1940];
    oli.NA = [99.6 170 251];
    oli.TMRCA = [1673 1710 1747];
end

% HCV Pybus mean parameters
x = oli.x(2);
y = oli.y(2);
NA = oli.NA(2);
NC = oli.NC(2);
r = oli.r(2);
TMRCA = oli.TMRCA(2);

% Mean stat of all the parameters
xm = mean(xest);
ym = mean(yest);
rm = mean(rest);
NAm = mean(NAest);
NCm = mean(NCest);
TMRCAm = mean(TMRCAest);

% Get Pybus population at all times - NA won't be consistent so recalculate
% as the means of estimates don't share the relations necessarily
NAfix = NC*exp(-r*(x-y));
Npyb = zeros(size(t));
for ii = 1:nTree
    tu = t(ii, :);
    I1 = tu >= x;
    I2 = tu > y & tu < x;
    I3 = tu <= y;
    Npyb(ii, :) = NC.*I1 + NAfix.*exp(-r.*(y - tu)).*I2 + NAfix.*I3;
end

% Get square error percent in population values
e = ((Nmean - Npyb)./Npyb).^2;

% Plot the estimates, their means and the pybus value
xtr = 1:nTree;
oxtr = ones(1, nTree);

figure;
subplot(6, 1, 1);
plot(xtr, xest, 1:nTree, x*oxtr, 1:nTree, xm*oxtr, 'k');
title(['Estimates of x with mean = [' num2str(xm) ' ' num2str(x) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('x');
subplot(6, 1, 2);
plot(xtr, yest, 1:nTree, y*oxtr, 1:nTree, ym*oxtr, 'k');
title(['Estimates of y with mean = [' num2str(ym) ' ' num2str(y) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('y');
subplot(6, 1, 3);
plot(xtr, rest, 1:nTree, r*oxtr, 1:nTree, rm*oxtr, 'k');
title(['Estimates of r with mean = [' num2str(rm) ' ' num2str(r) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('r');
subplot(6, 1, 4);
plot(xtr, NCest, 1:nTree, NC*oxtr, 1:nTree, NCm*oxtr, 'k');
title(['Estimates of NC with mean = [' num2str(NCm) ' ' num2str(NC) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('NC');
subplot(6, 1, 5);
plot(xtr, NAest, 1:nTree, NA*oxtr, 1:nTree, NAm*oxtr, 'k');
title(['Estimates of NA with mean = [' num2str(NAm) ' ' num2str(NA) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('NA');
subplot(6, 1, 6);
plot(xtr, TMRCAest, 1:nTree, TMRCA*oxtr, 1:nTree, TMRCAm*oxtr, 'k');
title(['Estimates of TMRCA with [mean pyb] = [' num2str(TMRCAm) ' ' num2str(TMRCA) ']']);
legend('est', 'pyb', 'mean', 'location', 'best');
xlabel('trees');
ylabel('TMRCA');

% Plot the space of N(t) estimates against Pybus2003


% Save batch data
save(['batchGarli_' num2str(nold(1))]);