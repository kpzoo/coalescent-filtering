% More tree testing code with the aim of drawing out the coalescent times
clear all
clc
close all

% Boolean to control plots
wantCoalTree = 0;

% Load data and get sequences
mu = 0.79*(10^-3);
name = 'result63.fasta';
seqs = fastaread(name);

% Take subset of sequences and get distances
num = 10;
seqs = seqs(1:num);
seqs = multialign(seqs);
%d = seqpdist(seqs, 'Alphabet', 'NT');
d = seqpdist(seqs, 'Alphabet', 'NT', 'Method', 'alignment-score');

% Get the tree and newick string
tr = seqlinkage(d, 'average', seqs);
tr = reroot(tr);
nu = getnewickstr(tr);

% Look at tree and get matrix and no. leaves, len = no. leaves and branches
phytreeviewer(tr);
[M, ID, D] = getmatrix(tr);
n = get(tr, 'NumLeaves');
len = 2*n - 1;

% Separate into branches and leaves
db = D(n+1:end);
dl = D(1:n);
idb = ID(n+1:end);
idl = ID(1:n);

% Get true branch lengths using the connectivity matrix and distances
[T, dcomp] = getCoalDists(M, n, D);

% Finally obtain full coalescent time vector
tcoal = [0 T];

if wantCoalTree
    % Reconstruct the tree from the distances by first taking relevant part of
    % M matrix out
    J = full(M);
    J = J(n+1:end, :);
    % Get children id for each branch
    idchild = zeros(n-1, 2);
    for i = 1:n-1
        idchild(i, :) = find(J(i, :));
    end
    
    % Reconstructed tree
    trec = phytree(idchild, T');
    phytreeviewer(trec);
end