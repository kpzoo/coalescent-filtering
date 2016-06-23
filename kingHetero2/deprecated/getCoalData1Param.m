% Function to generate coalescent data from a specific rate function

% Assumptions and Modifications
% - assumes sine function


function [twait, tcoal] = getCoalData1Param(fac, x, nData, w, A, lamtype)

% Define maximum rate L for sinusoid
%Lset = fac(1)*3*A*ones(size(fac));
Lset = fac*(x + A);

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
    % Calculate rate at current time
    switch(lamtype)
        case 1
            % Sinusoidal theta(t)
            lamt = fac(I)*(x*sin(w*t) + A);
        case 2
            % Sinusoidal N(t)
            lamt = fac(I)/(x*sin(w*t) + A);
    end
            
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