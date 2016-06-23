% Function to get the evolutionary distances in coalescent time format from
% the trees constructed vis seqlinkage
function [T, dcomp] = getCoalDists(M, n, D)

% Assumptions and modifications
% - calculates distance using both children and gives both for comparison
% - expects sparse M input as connectivity matrix of tree

% Get true branch lengths using the connectivity matrix and distances
T = -ones(1, n-1);
M = full(M);
len = 2*n - 1;
dcomp = -ones(2, n-1);

% First get distances for those branches that terminate with a leaf
for i = n+1:len
    % Get children ids that are leaves
    idz = find(M(i, 1:n));
    % Use sequence distance as branch length in this case
    if ~isempty(idz)
        T(i-n) = max(D(idz));
        dcomp(i-n, 1) = T(i-n);
        dcomp(i-n, 2) = min(D(idz));
    end
end

% Now get distances when the branch has only internal nodes as children
idnodist = find(T == -1) + n;
for i = 1:length(idnodist)
    % Children ids that are other branches ref to M total ids (2n-1)
    idbranch = find(M(idnodist(i), n+1:end)) + n;
    % Distance calculated based on children length and branch length
    dcomp(idnodist(i)-n, 1:2) = [D(idbranch(1)) + T(idbranch(1) - n) D(idbranch(2)) + T(idbranch(2) - n)];
    T(idnodist(i)-n) = max(dcomp(idnodist(i)-n, 1:2));
end



% Distance from subtracting branch distances in opposite direction then
    % adding the offset from previous node that is a child
    %T(idnodist(i)-n) = D(idbranch(1)) - D(idnodist(i)) + T(idbranch(1) - n);