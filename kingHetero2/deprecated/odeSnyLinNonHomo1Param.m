% Function to setup linear Snyder ODE set for use with ODE solvers

% Assumptions
% - modified for non-homogeneous simulations
% - only works with 1 RV and sine for now <-----------------------
% - y is a column vector, ts a vector just included for ode113
% - Q = [] means solving for random variable instead of DSPP
% - binfac is the appropriate binomial factor based on no. events

function dy = odeSnyLinNonHomo1Param(ts, y, Q, binfac, xset, w, lamtype)

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Rate function type
switch(lamtype)
    case 1
        % Sinusoidal theta(t)
        lamt = binfac*diag(xset*sin(w*ts) + max(xset));
    case 2
        % Sinusoidal N(t)
        lamt = binfac*diag(1./(xset*sin(w*ts) + max(xset)));
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