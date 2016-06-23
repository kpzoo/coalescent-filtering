% Function to perform Snyder ODE for exponential population function
function dy = odeSnyExp(ts, y, binfac, xsetMx)

% Assumptions
% - non-linear form with no DSPP

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for various functional forms
lamtdiag = binfac./(xsetMx(1, :).*exp(-xsetMx(2, :)*ts));

% Solve nonlinear differential equation set - RV filtering
nonLinDiag = lamtdiag*y;
dy = y'.*(-lamtdiag + nonLinDiag);

% Ensure output is column vector assuming input was
dy = dy';