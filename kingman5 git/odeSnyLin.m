% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP

function dy = odeSnyLin(ts, y, Q, lam)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Solve linear differential equation set
if isempty(Q)
    % RV filtering
    dy = y'*(-lam);
else
    % DSPP filtering
    dy = y'*(Q - lam);
end
% Ensure output is column vector assuming input was
dy = dy';