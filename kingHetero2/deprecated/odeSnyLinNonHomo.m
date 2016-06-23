% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - modified for non-homogeneous simulations
% - only works with 2 RVs for now <-----------------------
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyLinNonHomo(ts, y, Q, binfac, fntype, mi, xset)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Ensure correct no. of RVS
if length(mi) > 2
    error('ODE function can only handle up to 2 RVs');
end

% Inhomogeneous function type
switch(fntype)
    case 2
        % Exponential growth - let order be [x1(1) x2space], [[x1(1)
        % x2space]... etc
        lamt = zeros(1, prod(mi));
        % Iterate around variables to get order right
        id = 1;
        for i = 1:mi(1)
            nextid = id + mi(2);
            % Exponential function in partitions of param space
            lamt(id:nextid-1) = binfac*xset{1}(i)*exp(-xset{2}*ts);
            id = nextid;
        end
        % Convert to diagonal matrix of rates
        lamt = diag(lamt);
end

% Solve linear differential equation set
if isempty(Q)
    % RV filtering
    dy = y'*(-lamt);
else
    % DSPP filtering
    dy = y'*(Q - lamt);
end
% Ensure output is column vector assuming input was
dy = dy';