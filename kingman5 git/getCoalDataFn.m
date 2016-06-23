% Function to generate coalescent data from different functional forms of
% N(t) based on inputs

% Assumptions and Modifications
% - assumes either rejection sampling or time rescaling
% - includes logistic function and Pybus exponential

%% Main function to get coalescent data for various functions
function [twait, tcoal] = getCoalDataFn(fn)

% Get input data from structure
fnid = fn.id;
fac = fn.fac;
x = fn.param;
n = fn.nData + 1;

switch(fnid)
    case 1
        % Two parameter exponential - time rescaling
        % Simulate coalescent times for exponential using time rescaling
        tcoal = zeros(1, n);
        for i = 1:n-1
            % Get coalescent time for i+1 -> i lineages
            tcoal(i+1) = log(exp(x(2)*tcoal(i)) + x(2)*x(1)*exprnd(1)/fac(i))/x(2);
        end
        % Get waiting times
        twait = diff(tcoal);
        
    case 2
        % Four parameter sine - rejection
        [twait, tcoal] = getSinCoal(x, fac, fn.nData);
        
    case 3
        % One parameter constant - direct from exponential waiting times
        rates = fac/x(1);
        twait = exprnd(1./rates);
        tcoal = cumsum([0 twait]);
    case 4
        % Four parameter logistic - rejection
        [twait, tcoal] = getLogCoal(x, fac, fn.nData, fnid);
    case 5
        % Four parameter Pybus exponential - rejection
        [twait, tcoal] = getLogCoal(x, fac, fn.nData, fnid);
end





%% Subfunction for sinsusoidal

% Function to generate N(t) = x1sin(x2t + x3) + x4 data by rejection 
function [twait, tcoal] = getSinCoal(x, fac, nData)

% Ensure x of correct length
if length(x) ~= 4
    error('Expect 4 parameters');
end

% Define maximum rate L for sinusoid N(t)
%Lset = fac/(x(4) - x(1));
Lset = (fac(1)/(x(4) - x(1)))*ones(size(fac));

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
    lamt = fac(I)/(x(1)*sin(x(2)*t + x(3)) + x(4));
            
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


%% Subfunction for logistic

% Function to generate N(t) = x1(1 + exp(-x2x3))/(1 + exp(-x2(x3 - t)) + x4 data by rejection 
function [twait, tcoal] = getLogCoal(x, fac, nData, fnid)

% Ensure x of correct length
if length(x) ~= 4
    error('Expect 4 parameters');
end

% Define maximum rate L for sinusoid N(t)
%Lset = fac/x(4);
Lset = (fac(1)/x(4))*ones(size(fac));

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
    
    % Calculate rate at current time for logistic or Pybus N(t)
    switch(fnid)
        case 4
            % Logistic
            Nt = x(1)*(1 + exp(-x(2)*x(3)))/(1 + exp(-x(2)*(x(3) - t))) + x(4);
        case 5
            % Pybus piecewise exponential using indicators
            I1 = t <= x(3);
            I2 = t > x(3) && t < x(4);
            I3 = t >= x(4);
            Nt = x(1)*I1 + x(1)*exp(-x(2)*(t - x(3)))*I2 + ...
                x(1)*exp(-x(2)*(x(4) - x(3)))*I3;
    end
    lamt = fac(I)/(Nt);
            
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