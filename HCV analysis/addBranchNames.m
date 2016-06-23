% Basic function to take a newick string and add branch names
function newStr = addBranchNames(newIn, newOutName, writeYes)

% Assumptions and modifications
% - takes input as tree or newick string for newIn

% Check input
if ~ischar(newOutName)
    error('Did not specify output string for writing to file');
end

% Read tree if given string else input is tree
if ischar(newIn)
    tree = phytreeread(newIn);
else
    tree = newIn;
end

% Write branch names
if writeYes
    fullname = [newOutName '.tree'];
    phytreewrite(fullname, tree, 'Branchnames', 'true');
end

% Get newick string for tree with branch names
newStr = getnewickstr(tree, 'Branchnames', 'true');