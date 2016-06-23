% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - modified for non-homogeneous simulations of N = N0exp(-rt) but with
% both parameters to be estimated <----------------------
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyExp2Param(ts, y, Q, ldiagTno, xset, mi)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end
if length(mi) ~= 2
    error('Function expects 2 random variables');
end

% Construct time dependent rate matrix diagonals by noting that lamnoTdiag
% must have mi(2) repetitions of mi(1) values
m = prod(mi);
lamdiag = zeros(1, m);
for kk = 1:mi(2)
    % Index the start and end points based on mi(2)
    id1 = 1 + (kk - 1)*mi(1);
    id2 = kk*mi(1);
    x2 = xset{2}(kk);
    lamdiag(1, id1:id2) = ldiagTno(1, id1:id2)*exp(x2*ts);
end

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