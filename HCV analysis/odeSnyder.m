% Function to setup linear Snyder ODE set for use with ODE solvers for RVs

% Assumptions
% - removed DSPP allowance
% - to improve speed removed y'*diag(lamtdiag) with y'.*lamtdiag
% - generalised for non-homogeneous simulations of N(t)
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyder(ts, y, binfac, fn, xsetMx)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for various functional forms
[~, lamtdiag] = getTimeVaryingN(fn, ts, binfac, xsetMx);

% Solve linear differential equation set - RV filtering
dy = y'.*(-lamtdiag);

% Ensure output is column vector assuming input was
dy = dy';