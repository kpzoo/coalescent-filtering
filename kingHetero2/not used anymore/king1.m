%% Basic Kingman coalescent code with process statistics

% Assumptions
% - constant population size
% - isochronous sampling
% - neutral mutations

clc
close all
clear all

% Set sample parameters
gam = 0.57721;
n = 5;
kset = n:-1:2;
klen = length(kset);
nRuns = 1;
tkset = zeros(klen, nRuns);
rateset = kset.*(kset-1)/2;

%% Simulate the coalescent event times

% Draw the coalescent times based on the rates
for i = 1:klen
    mu = 1/rateset(i);
    tkset(i, 1:nRuns) = exprnd(mu, [1 nRuns]);
end

% Get MRCA and total time samples
tmrcaset = sum(tkset);
ttotset = sum(tkset'*(kset'), 2);
t2set = tkset(klen, :);

% Calculate time statistics - empirical and theoretical
tmrca.avg = mean(tmrcaset);
tmrca.avgtheo = 2*(1 - 1/n);
tmrca.var = var(tmrcaset);
tmrca.vartheo = 8*sum(1./(kset.*kset)) - tmrca.avgtheo^2;

ttot.avg = mean(ttotset);
ttot.avgtheo = 2*(gam + log(n));
ttot.var = var(ttotset);
ttot.vartheo = 4*sum(1./((kset-1).*(kset-1)));

t2.avg = mean(t2set);
t2.avgratio = t2.avg/tmrca.avg;
t2.var = var(t2set);
t2.varratio = t2.var/tmrca.var;

%% Markov death process for the coalescent no. of equivalence classes |R|

rr = 2:n;
% Define the death rates - using the function f(|R|) from Kingman with
% state 1 as absorbing state
Qlow = rr.*(rr - 1)/2;                                                                            
Q = diag(Qlow, -1) + diag([0 -Qlow]);


%% Generate a random tree with the coalescent times as distances

% Define storage cells
Tree = cell(1, nRuns);
Children = cell(1, nRuns);
Labels = cell(1, nRuns);

% Branch length cell which is used for mutation count and mutation rate, mu
C = cell(1, nRuns);
mu = 1;
nMutations = zeros(nRuns, n-1);

% Calculate a tree with coalescent distances for each run
for i = 1:nRuns
    % Get the metric tree and the branch to root distances C by summing the
    % inter-coalescent times
    C{i} = cumsum(tkset(:, i));
    [Children{i}, Labels{i}, Tree{i}] = getKingmanTree(n, C{i});
    
    % Get the number of mutations for each branch of each tree
    lbranch = C{i}';
    rbranch = mu*lbranch;
    % Rate for Poisson mutations based on branch lengths
    nMutations(i, :) = poissrnd(rbranch);
end

%% Analyse a tree of choice

% Get properties of the mth tree
m = 1;
tree = Tree{m};
% Gives sparse connection matrix of tree
[M, ID, Dists] = getmatrix(tree);
phytreeviewer(tree);

