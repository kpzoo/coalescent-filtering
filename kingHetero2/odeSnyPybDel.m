% Function to perform Snyder ODE for Pybus population function in which the
% x4 parameter is the period of exponential so that x3+x4 is the end time
function dy = odeSnyPybDel(ts, y, binfac, xsetMx, delRV)

% Assumptions
% - non-linear form with no DSPP

% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for various functional forms
cond1 = delRV.x3;
cond2 = cond1 + delRV.x4;
I1 = ts <= cond1;
I2 = ts > cond1 & ts < cond2;
I3 = ts >= cond2;
Ntemp = xsetMx(1, :).*I1 + xsetMx(1, :).*exp(-xsetMx(2, :).*(ts - delRV.x3)).*I2 + ...
    xsetMx(1, :).*exp(-xsetMx(2, :).*delRV.x4).*I3;
lamtdiag = binfac./Ntemp;

% Solve nonlinear differential equation set - RV filtering
nonLinDiag = lamtdiag*y;
dy = y'.*(-lamtdiag + nonLinDiag);

% Ensure output is column vector assuming input was
dy = dy';