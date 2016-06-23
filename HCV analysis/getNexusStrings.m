% Extract newick strings from a nexus file with several trees
function [nTree, treeStr] = getNexusStrings(name)

% Assumptions and modifications
% - assumes the text file is from the figtree file produced by garli
% - assumes structure tree rep1 = [stuff] (newick)

% Read data from the input text file into a cell of strings
 A = importdata(name);
 numStr = length(A);
 j = 0;
 treeStr = cell(1, 1);
 
 % Loop across data and find the trees and take the newick form
 for i = 1:numStr
     % Find all the lines with trees in them
     stind = regexp(A{i}, 'tree rep', 'once');
     % Extract the newick
     if ~isempty(stind)
         dataStr = A{i};
         % Starting index of newick part
         st = regexp(dataStr, '(');
         % Store tree and update cell
         j = j + 1;
         % Assumes rest of line from start index is newick string
         treeStr{j} = dataStr(st(1):length(dataStr));
     end
 end
 
% Number of trees in file
nTree = length(treeStr);