% Function to generate N(t) = x1sin(x2t + x3) + A data by rejection 
function [twait, tcoal] = getSinCoal(x, fac)

% Ensure x of correct length
if length(x) ~= 3
    error('Expect 3 parameters');
end

% Define maximum rate L for sinusoid N(t)
%Lset = fac/(A - x);
Lset = (fac(1)/(A - x(1)))*ones(size(fac));

% Duplicate smallest L value for nData+1 case assuming falling binomial fac
Lset = [Lset Lset(end)];

% Simulate process via thinning algorithm assuming lam(t) <= L for t <= T
I = 1;
t = 0;
tcoal = zeros(1, nData+1);
while(I <= nData) 
    % Generate a Poisson homogeneous interval
    U = rand;
    t = t -log(U)/Lset(I);
    % Calculate rate at current time for sinusoidal N(t)
    lamt = fac(I)/(x(1)*sin(x(2)*t) + x(3));
            
    % Rejection sample by calculating rate
    U = rand;
    if U <= lamt/Lset(I)
        % An event has occurred so take data, ensure tcoal(1) = 0
        I = I + 1;
        tcoal(I) = t;
    end
    %disp(['Finished event ' num2str(I)]);
end

% Obtain interval data
twait = diff(tcoal);