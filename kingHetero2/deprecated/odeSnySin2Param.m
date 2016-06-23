% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - modified for non-homogeneous simulations of N(t) = A + x1sin(x2t)
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnySin2Param(ts, y, Q, binfac, xset, mi, A)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end
if length(mi) ~= 2
    error('Function expects 2 random variables');
end

% Construct time dependent rate matrix diagonals using kronecker products
% of x1sin(x2t)
x1 = xset{1};
x2 = sin(xset{2}*ts);
Nt = kron(x1, x2) + A;
lamdiag = binfac./Nt;

% Rate matrix
lamt = diag(lamdiag);

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