%%%%% ERRORR with how handling Lset


% Function to generate coalescent data from a specific rate function

function [twait, tcoal] = getCoalescentData(fntype, fac, x, nData)

% Assumptions and Modifications
% - the Lset which controls the rejection must be set for each function



% List functional forms of N(t) and expected no. RVs and check input
fnset = {'const', 'expd'};
fnRV = [1 2];
if length(x) ~= fnRV(fntype)
    assignin('base', 'noInp', length(x));
    assignin('base', 'noExp', fnRV(fntype));
    error('Mismatch between function and expected no. of RVs');
end

% Define maximum rate L based on function and max coalescent binomial fac
switch(fntype)
    case 1
        % Constant value - assumes working in lam = 1/N
        Lset = fac*x;
    case 2
        % Exponential decay
        Lset = fac*x(1);
end
% Duplicate smallest L value for nData+1 case assuming falling binomial fac
Lset = [Lset Lset(end)];

% Simulate process via thinning algorithm assuming lam(t) <= L for t <= T
I = 1;
t = 0;
tcoal = zeros(1, nData+1);
while(I <= nData)
    % Calculate rate at current time and time vary Lset is needed
    switch(fntype)
        case 1
            % Constant
            lamt = fac(I)*x;
            Lset(I) = fac(I)*x;
        case 2
            % Exponential growth
            lamt = fac(I)*x(1)*exp(-x(2)*t);
            Lset(I) = fac(I)*x(1);
    end
    
    % Generate a Poisson homogeneous interval
    U = rand;
    t = t -log(U)/Lset(I);
    
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