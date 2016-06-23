% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - modified for non-homogeneous simulations of N = N0exp(-xt)
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyExp1Param(ts, y, Q, binfac, xset, N0)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Rate function type
lamt = diag((binfac/N0)*exp(xset*ts));

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