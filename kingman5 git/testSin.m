% Assumes batch data available and in workspace

% Draw out parameters
y = cell(1, 1);
x = zeros(1, numRV);
for i = 1:numRV
    y{i} = paramhat(i, :);
    x(i) = param(i);
end

% Parse the batch data based on function type
switch(fn.id)
    case 2
        % Sinusoidal function
        
        % Get cycle and make t sensible at 100 cycles with 50 points per
        period = 2*pi/x(2);
        maxt = 100*period;
        dt = period/50;
        
        % Reconstruct N(t) and it's estimates
        t = 0:dt:maxt;
        t = t';
        lent = length(t);
        Nhat = zeros(lent, M);
        N = x(1)*sin(x(2)*t + x(3)) + x(4);
        
        % Get the N(t) estimated waveforms
        for i = 1:lent
            for j = 1:M
                Nhat(i, j) = y{1}(j)*sin(y{2}(j)*t(i) + y{3}(j)) + y{4}(j);
            end
        end
        
    case 4
        % Logistic function
        
        % Make t sensible at twice the value of t50 = x3
        dt = 1;
        maxt = x(3)*2;
        t = 0:dt:maxt;
        t = t';
        
        % Reconstruct N(t) and it's estimates
        N = x(1)*(1 + exp(-x(2)*x(3)))./(1 + exp(-x(2)*(x(3) - t))) + x(4);
        lent = length(t);
        Nhat = zeros(lent, M);
        
        % Get the N(t) estimated waveforms
        for i = 1:lent
            for j = 1:M
                Nhat(i, j) = y{1}(j)*(1 + exp(-y{2}(j)*y{3}(j)))/(1 + exp(-y{2}(j)*(y{3}(j) - t(i))))...
                    + y{4}(j);
            end
        end
        
    case 5
        % Piecewise exponential
        
        % Make t sensible based on x = x(3) and y = x(4)
        maxt = x(4) + x(3);
        dt = maxt/10000;
        t = 0:dt:maxt;
        t = t';
        
        % Reconstruct N(t) and it's estimates
        I1 = t <= x(3);
        I2 = t > x(3) & t < x(4);
        I3 = t >= x(4);
        N = x(1).*I1 + x(1)*exp(-x(2)*(t - x(3))).*I2 + ...
            x(1)*exp(-x(2)*(x(4) - x(3))).*I3;
        lent = length(t);
        Nhat = zeros(lent, M);
        
        
        % Get the N(t) estimated waveforms
        for i = 1:lent
            for j = 1:M
                I1 = t(i) <= y{3}(j);
                I2 = t(i) > y{3}(j) & t(i) < y{4}(j);
                I3 = t(i) >= y{4}(j);
                Nhat(i, j) = y{1}(j)*I1 + y{1}(j)*exp(-y{2}(j)*(t(i) - y{3}(j)))*I2 + ...
                    y{1}(j)*exp(-y{2}(j)*(y{4}(j) - y{3}(j)))*I3;
            end
        end 
end

% Get correlation
pmcc = zeros(M, 1);
for i = 1:M
    pmcc(i) = corr2(N, Nhat(:, i));
end

% Percentage square errors
mse = zeros(M, 1);
for i = 1:M
    mse(i) = mean(100*((Nhat(:, i) - N).^2)./(N.^2));
end

% Plot stats
Mmse = mean(mse)*ones(size(mse));
Smse = mean(mse)*ones(size(mse));
mseUB = (Mmse + 2*Smse);
mseLB = (max(0, Mmse - 2*Smse));
figure;
plot(1:M, mse, 1:M, Mmse, 'k', 1:M, mseLB, 'r', 1:M, mseUB, 'r');

% Get histogram
figure;
hist(mse);
xlabel('% relative MSE');
ylabel('frequency');
title(['Performance across runs, E[mse] = ' num2str(mean(mse))]);