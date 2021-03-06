% Function to get optimal estimate of N(t) and its bounds based on the
% final posterior distribution and going forward in time
function [Nf, Nbnd, qcum] = getPopForwardQuan(t, q, fn, otherPms)

% Assumptions and modifications
% - includes quantile calculation instead of just std as better maybe
% - includes case 7 
% - removed a loop to make more efficient
% - works with time going forwards (not usual backward coalescent)
% - q is the final posterior and t is vector of times

% Extract data of interest
pm = fn.xsetMx;
fnid = fn.id;
lent = length(t);
Nf = zeros(1, lent);
Nbnd = zeros(2, lent);
%Nquan = zeros(2, lent);

% Check that posterior and inputs same length
lenq = length(q);
if size(pm, 2) ~= lenq
    error('Inputs inconsistent in dimension');
end

% Quantiles to cover 95% confidence
qcum = cumsum(q);
id2p5 = find(qcum <= 0.025, 1, 'last');
id97p5 = find(qcum >= 0.975, 1, 'first');
id50s(1) = find(qcum <= 0.5, 1, 'last');
id50s(2) = find(qcum >= 0.5, 1, 'first');
did = abs(0.5 - qcum(id50s));
id50 = id50s(did == min(did));

% Loop across times and input parameters for model
for i = 1:lent
     % Pick function and calculate population space at each time
     switch(fnid)
         case 5
             % Piecewise exponential in forward time Pybus2003
             tpres = otherPms;
             x = tpres - pm(3, :);
             y = tpres - pm(4, :);
             NC = pm(1, :);
             r = pm(2, :);
             NA = NC.*exp(-r.*(pm(4, :)-pm(3, :)));
             
             % Indicator functions and population
             tz = t(i);
             I1 = tz >= x;
             I2 = tz > y & tz < x;
             I3 = tz <= y;
             Nq = NC.*I1 + NA.*exp(-r.*(y - tz)).*I2 + NA.*I3;
             
         case 6
             % Piecewise exponential in forward time reparameterised
             tpres = otherPms;
             x = tpres - pm(3, :);
             % The y is the diff param + backward time x then use tpres
             y = tpres - (pm(4, :) + pm(3, :));
             NC = pm(1, :);
             r = pm(2, :);
             NA = NC.*exp(-r.*pm(4, :));
             
             % Indicator functions and population
             tz = t(i);
             I1 = tz >= x;
             I2 = tz > y & tz < x;
             I3 = tz <= y;
             Nq = NC.*I1 + NA.*exp(-r.*(y - tz)).*I2 + NA.*I3;
             
         case 7
             % Piecewise exponential in log form with x and y-x
             tpres = otherPms;
             x = tpres - pm(3, :);
             % The y is the diff param + backward time x then use tpres
             y = tpres - (pm(4, :) + pm(3, :));
             NC = exp(pm(1, :));
             r = pm(2, :);
             NA = NC.*exp(-r.*pm(4, :));
             
             % Indicator functions and population
             tz = t(i);
             I1 = tz >= x;
             I2 = tz > y & tz < x;
             I3 = tz <= y;
             Nq = NC.*I1 + NA.*exp(-r.*(y - tz)).*I2 + NA.*I3;
     end
      
     % Probability statistics at each time
     Nf(i) = q*Nq';
     Nqstd = sqrt(q*((Nq.^2)') - Nf(i)^2);
     Nbnd(:, i) = [Nf(i) - 2*Nqstd, Nf(i) + 2*Nqstd];
     Nbnd(Nbnd < 0) = 0;
     %Nf(i) = Nq(id50);
     
     % Percentiles for 95% cover
     %Nbnd(:, i) = [Nq(id2p5) Nq(id97p5)];
     %Neg = Nq(id50);
     %Nbnd(:, i) = [Nq(id2p5) Nq(id97p5)] - Neg;
     %Nbnd(:, i) = Nbnd(:, i) + Nf(i);
end