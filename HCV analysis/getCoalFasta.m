% Function to read fasta files, create and root a tree and then get
% coalescent times
function [tcoal, twait, trees, n, evDist, dcomp] = getCoalFasta(mu, name)

% Assumptions and modifications
% - properly reconstructs the timing tree from the coalescent times
% - NOT SURE IF CONVERSION OF BRANCHES TO TIMES IS CORRECT
% - the name file must be a .fasta and in the current folder

% Read the fasta data sets
seqs = fastaread(name);
%seqs = multialign(seqs);

% Obtain pairwise genetic distances and create a tree
%distSeq = seqpdist(seqs, 'Alphabet', 'NT', 'Method', 'Hasegawa', 'Indels','complete-delete');
distSeq = seqpdist(seqs, 'Alphabet', 'NT');
treeSeq = seqlinkage(distSeq, 'weighted', seqs);

% Reroot the trees via the midpoint method and get newick string
%treeSeq = reroot(treeSeq);
newickSeq = getnewickstr(treeSeq);
trees.treeSeq = treeSeq;
trees.newickSeq = newickSeq;

% Get the branch to leaf distances and the no. lineages
[M, ~, D] = getmatrix(treeSeq);
n = get(treeSeq, 'NumLeaves');

% Get true branch lengths using the connectivity matrix and distances
[T, dcomp] = getCoalDists(M, n, D);

% Finally obtain full coalescent time vector - converting distance to times
evDist = [0 T];
evDist = sort(evDist);
tcoal = evDist/mu; %% <-----------------------------
twait = diff(tcoal);

% Reconstruct the tree from the distances by first taking relevant part of
% M matrix out
J = full(M);
J = J(n+1:end, :);
% Get children id for each branch
idchild = zeros(n-1, 2);
for i = 1:n-1
    idchild(i, :) = find(J(i, :));
end
tree = phytree(idchild, T');

% Output data
newick = getnewickstr(tree);
trees.treeCoal = tree;
trees.newickCoal = newick;
trees.idchild = idchild;