% Function to calculate time dependent rate matrix diagonals
function [Nt, lamtdiag] = getTimeVaryingN(fn, t, binfac, xGen)

% Assumptions and modifications
% - includes logistic as case 4 and Pybus exponential as case 5
% - allows calculation of true value when specified with true parameters
% - the fn set must be consistent with the id
% - time input must be a scalar
% - removed trueVal boolean and replaced with xGen which can be xsetMx or x

% % Check time input is scalar
% if ~all(size(t) == 1)
%     error('t needs to be a scalar');
% end
    
% Extract inputs
id = fn.id;
numRV = fn.numRV;
m = fn.m;

% Get variables in form for direction function input
y = cell(1, numRV);
for i = 1:numRV
        % Want space of Nt so use xsetMx or column form of x if want true
        % value instead of estimates
        y{i} = xGen(i, :);
end

% Choose functional form and calculate population function
switch(id)
    case 1
        % 2 parameter exponential
        Nt = y{1}.*exp(-y{2}*t);
    case 2
        % 4 parameter sinusoidal N(t)
        Nt = y{1}.*sin(y{2}*t + y{3}) + y{4};
    case 3
        % Constant N(t)
        Nt = y{1};  
    case 4
        % Logistic N(t)
        Nt = y{1}.*(1 + exp(-y{2}.*y{3}))./(1 + exp(-y{2}.*(y{3} - t))) + y{4};
    case 5
        % Pybus piecewise exponential N(t)
        I1 = t <= y{3};
        I2 = t > y{3} & t < y{4};
        I3 = t >= y{4};
        Nt = y{1}.*I1 + y{1}.*exp(-y{2}.*(t - y{3})).*I2 + ...
            y{1}.*exp(-y{2}.*(y{4} - y{3})).*I3;
end

% Get the rate matrix diagonal from coalescent form and check dimension -
% note that when trueVal == 1 it is no longer a diagonal
lamtdiag = binfac./Nt;
% actualDim = length(lamtdiag);
% if actualDim ~= m && ~trueVal
%     assignin('base', 'actualDim', actualDim);
%     assignin('base', 'expDim', m);
%     error('Incorrect diagonal assembled');
% end