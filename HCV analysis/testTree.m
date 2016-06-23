% Code to create a tree from the data of Pybus 2003 and then collect the
% coalescent times for use with the Snyder filter
clear all
clc
close all

% Read the fasta data sets
seqs1 = fastaread('result63.fasta');
seqs2 = fastaread('result68.fasta');

% Obtain pairwise distances and create a tree
distances1 = seqpdist(seqs1, 'Method', 'Jukes-Cantor');
tree1 = seqlinkage(distances1, 'UPGMA', seqs1);
distances2 = seqpdist(seqs2, 'Method', 'Jukes-Cantor');
tree2 = seqlinkage(distances2, 'UPGMA', seqs2);

% Reroot the trees via the midpoint method
tree1 = reroot(tree1);
tree2 = reroot(tree2);

% Get the newick string of each tree
str1 = getnewickstr(tree1);
str2 = getnewickstr(tree2);

% View the tree in evolutionary distance for both data sets
h1 = plot(tree1, 'orient', 'top');
ylabel('Evolutionary distance')
set(h1.terminalNodeLabels, 'Rotation', 65)
h2 = plot(tree2, 'orient', 'top');
ylabel('Evolutionary distance')
set(h2.terminalNodeLabels, 'Rotation', 65)

% Get the distances from the tree and convert to times
[M1, ID1, Dists1] = getmatrix(tree1);
[M2, ID2, Dists2] = getmatrix(tree2);
mu = 0.79*(10^-3);

% The distances extracted are for branches and leaves - only want branches
numLeaves1 = get(tree1, 'NumLeaves');
numLeaves2 = get(tree2, 'NumLeaves');
dbranch1 = Dists1(numLeaves1+1:end);
dbranch2 = Dists2(numLeaves2+1:end);

% Convert branch distances to times
t1 = dbranch1/mu;
t2 = dbranch2/mu;

% The times are with respect to the leaves but out of order
t1 = sort(t1);
t2 = sort(t2);

% See if can reconstruct the tree sensibly with these times
[children, labels, tree] = getKingmanTree(numLeaves1, t1);
phytreeviewer(tree);
[children, labels, tree] = getKingmanTree(numLeaves2, t2);
phytreeviewer(tree);