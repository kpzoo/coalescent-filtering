% Function to create a Kingman coalescent tree

% Assumptions:
% - exchangability so labels don't matter except to name parents
% - the coalescing is random and the parent given a new id

function [children, labels, tree] = getKingmanTree(n, C)

% Inputs:
% n - number of lineages
% C - distances from branches to leaves (column vector)

% Outputs:
% children - the list of coalesced lineages by numbers from the labels
% tree - the phylogenetic tree object

% Check n and C are appropriately sized and if ultrametric
nBranches = n - 1;
if ~isempty(C)
    ultra = 0;
    if length(C) ~= nBranches
        error('Incorrect tree distance number');
    end
else
    % Ultrametric tree chosen
    ultra = 1;
end

% Every generation randomly pick 2 lineages to coalesce and store there
% identifiers while also updating the labels by +1 each generation
children = zeros(nBranches, 2);
labels = cell(1, n);
labels{1} = 1:n;
for i = 1:nBranches
    % Randomly choose 2 nodes to coalesce
    children(i, :) = randsample(labels{i}, 2);
    % Update labels with new id for coalesced parent
    remLabels = setdiff(labels{i}, children(i, :));
    newLabel = n + i;
    labels{i+1} = [remLabels newLabel];
end

% Create the phylogenetic tree object
if ultra
    tree = phytree(children);
else
    % Include distances based on coalescent process
    tree = phytree(children, C);
end