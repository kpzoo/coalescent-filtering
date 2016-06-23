% Function to get optimal estimate of N(t) and its bounds based on the
% final posterior distribution and going backward in time like coalescent
function [Nf, Nbnd] = getPopBackward(t, q, fn)

% Assumptions and modifications
% - modified to ensure no negative square bounds
% - removed a loop to make more efficient
% - modification of getPopForward that works in coalescent backward time
% - q is the final posterior and t is vector of times

% Extract data of interest
pm = fn.xsetMx;
lent = length(t);
Nf = zeros(1, lent);
Nbnd = zeros(2, lent);

% Check that posterior and inputs same length
lenq = length(q);
if size(pm, 2) ~= lenq
    error('Inputs inconsistent in dimension');
end

% Loop across times and input parameters for model
for i = 1:lent     
     % Select equation based on functional form
     switch(fn.id)
         case 1
             % 2 parameter exponential
             Nq = pm(1, :).*exp(-pm(2, :)*t(i));
         case 2
             % 4 parameter sinusoidal N(t)
             Nq = pm(1, :).*sin(pm(2, :)*t(i) + pm(3, :)) + pm(4, :);
         case 3
             % Constant N(t)
             Nq = pm(1, :);
         case 4
             % Logistic N(t)
             Nq = pm(1, :).*(1 + exp(-pm(2, :).*pm(3, :)))./...
                 (1 + exp(-pm(2, :).*(pm(3, :) - t(i)))) + pm(4, :);
         case 5
             % Pybus piecewise exponential N(t)
             I1 = t(i) <= pm(3, :);
             I2 = t(i) > pm(3, :) & t(i) < pm(4, :);
             I3 = t(i) >= pm(4, :);
             Nq = pm(1, :).*I1 + pm(1, :).*exp(-pm(2, :).*(t(i) - pm(3, :))).*I2 + ...
                 pm(1, :).*exp(-pm(2, :).*(pm(4, :) - pm(3, :))).*I3;
         case 6
             % Pybus piecewise exponential N(t) reparametrised to x and y-x
             I1 = t(i) <= pm(3, :);
             I2 = t(i) > pm(3, :) & t(i) < pm(4, :) + pm(3, :);
             I3 = t(i) >= pm(4, :) + pm(3, :);
             Nq = pm(1, :).*I1 + pm(1, :).*exp(-pm(2, :).*(t(i) - pm(3, :))).*I2 + ...
                 pm(1, :).*exp(-pm(2, :).*(pm(4, :))).*I3;
     end
     
    % Probability statistics at each time
    Nf(i) = q*Nq';
    Nqstd = sqrt(max(0, q*((Nq.^2)') - Nf(i)^2)); %%%<--- max(0, ) to ensure no negatives
    Nbnd(:, i) = [Nf(i) - 2*Nqstd, Nf(i) + 2*Nqstd];  
end