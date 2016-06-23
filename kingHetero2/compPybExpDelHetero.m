% Script to compare batch runs on 2 parameter pybus and exp N(t) that are
% heterochronously sampled: con-exp-con parametrised as in HCV

% Assumptions and Modifications
% -

clearvars
clc
close all

%% Set loop parameters and run simulation and filtering for exponential

% Initialise total number of samples and use factors to get nSampTimes
n = 200;
zz = 1:n;
nSampTimes = zz(~(rem(n, zz)));
nSampTimes = nSampTimes(1:end-1); % don't want nSampTimes = n case
sampSize = n./nSampTimes;

% Define no. runs and repetitions and maxSamp time and filter dimension
nReps = 5000;
nRuns = length(nSampTimes);
numRV = 2;
mi = 30*ones(1, numRV);

% Variables to store outputs - each row is a different sampling scheme
ss = zeros(nRuns, nReps);
x1 = ss; x1hat = ss; x2 = ss; x2hat = ss; x1a = ss; x1ahat = ss; x2a = ss; x2ahat = ss;
Nz = cell(nRuns, nReps); Nz1 = Nz; Nzhat = Nz; Nzhat1 = Nz;

% Main loop to repeat many iterations of the same sampling and then also
% loop across different sampling schemes
for ii = 1:nRuns
    for jj = 1:nReps
        % Main simulation of coalescent with filtering
        [est, est1, tz, Nz{ii, jj}, Nz1{ii, jj}, Nzhat{ii, jj}, Nzhat1{ii, jj}, extra] = pybexpDelFn(n, nSampTimes(ii), mi);
        % Parameters of interest
        x1(ii, jj) = est.x(1); x2(ii, jj) = est.x(2);
        x1a(ii, jj) = est1.x(1); x2a(ii, jj) = est1.x(2);
        x1hat(ii, jj) = est.xhat(1); x2hat(ii, jj) = est.xhat(2);
        x1ahat(ii, jj) = est1.xhat(1); x2ahat(ii, jj) = est1.xhat(2);
    end
    disp(['Finished sampling scheme ' num2str(ii) ' of ' num2str(nRuns)]);
end

% Collect data into more manageable and extendable format
err{1} = (x1 - x1hat)./x1; err{2} = (x2 - x2hat)./x2;
err1{1} = (x1a - x1ahat)./x1a; err1{2} = (x2a - x2ahat)./x2a;
meanE = cell(1, numRV); meanE1 = cell(1, numRV);
mseE = meanE; mseE1 = meanE1;
for i = 1:numRV
    esq = err{i}.^2;
    esq1 = err1{i}.^2;
    meanE{i} = mean(err{i}, 2);
    mseE{i} = mean(esq, 2);
    meanE1{i} = mean(err1{i}, 2);
    mseE1{i} = mean(esq1, 2);
end


% Performs simple ranking procedure - assign value nReps to worst MSE or
% mean error and 1 to lowest
MSErnks = zeros(numRV, nRuns); MSErnks1 = zeros(numRV, nRuns);
Meanrnks = MSErnks; Meanrnks1 = MSErnks1;
for i = 1:numRV
    MSErnks(i, :) = tiedrank(mseE{i});
    Meanrnks(i, :) = tiedrank(abs(meanE{i}));
    MSErnks1(i, :) = tiedrank(mseE1{i});
    Meanrnks1(i, :) = tiedrank(abs(meanE1{i}));
end

% Then get overall rank based on summing ranks of each scheme for each xi
schemeMSE = tiedrank(sum(MSErnks));
schemeMean = tiedrank(sum(Meanrnks));
schemeMSE1 = tiedrank(sum(MSErnks1));
schemeMean1 = tiedrank(sum(Meanrnks1));

% Convert demographic data to more usable form with MSE across a trajectory
% and mean absolute error
eN = zeros(nRuns, nReps); eN1 = eN;
eNsq = eN; eNsq1 = eN;
for i = 1:nRuns
    for j = 1:nReps
        e = Nzhat{i, j} - Nz{i, j};
        e1 = Nzhat1{i, j} - Nz1{i, j};
        eN(i, j) = trapz(tz, e)/range(tz);
        eN1(i, j) = trapz(tz, e1)/range(tz);
        eNsq(i, j) = trapz(tz, e.^2)/range(tz);
        eNsq1(i, j) = trapz(tz, e1.^2)/range(tz);
    end
end

% Check for fixed true trajectory
if all(Nz{1, 1} == Nz{nRuns, nReps}) && all(Nz1{1, 1} == Nz1{nRuns, nReps})
    comparison = 1;
    Nz = Nz{1, 1};
    Nz1 = Nz1{1, 1};
else
    comparison = 0;
end

% Control the x axis labelling and form
dom = 1:nRuns;
xlab = cell(1, nRuns);
for i = dom
    xlab{i} = num2str(nSampTimes(i));
end
xlab2 = cell(1, nRuns);
for i = dom
    xlab2{i} = ['K = ' num2str(nSampTimes(i))];
end


% Save the data
save(['bigpyb' num2str(n)]);


%% Plotting for the population functions

% Compare probability distributions from histograms for population MSE
figure;
for i = 1:nRuns
    subplot(ceil(nRuns/2), 2, i);
    h = histogram(eN(i, :), 20);
    h.Normalization = 'probability';
    hold all
    h = histogram(eN1(i, :), 20);
    h.Normalization = 'probability';
    hold off
    grid;
    xlabel(['<N - Nhat>']);
    ylabel('probability');
    legend('con-exp-con', 'exponential', 'location', 'best');
    title(['Population mean error for ' num2str(nSampTimes(i)) ' sample points']);
end

figure;
for i = 1:nRuns
    subplot(ceil(nRuns/2), 2, i);
    h = histogram(eNsq(i, :), 20);
    h.Normalization = 'probability';
    hold all
    h = histogram(eNsq1(i, :), 20);
    h.Normalization = 'probability';
    hold off
    grid;
    xlabel(['<N - Nhat>']);
    ylabel('probability');
    legend('con-exp-con', 'exponential', 'location', 'best');
    title(['Population MSE for ' num2str(nSampTimes(i)) ' sample points']);
end

% If only one trajectory in all of Nz then plot this and errors
if comparison
    
    % Take a subset of possible sample schemes and plot bounds
    nPlots = 4;
    iPlots = round(linspace(1, nRuns, nPlots));
    ax = zeros(1, nPlots); % to link axes
    % Con-exp-con function - bounds and estimate
    figure;
    Nreal = Nz;
    % Get standard deviation bounds
    sig = zeros(nPlots, length(tz));
    mu = sig;
    for i = 1:nPlots
        ax(i) = subplot(ceil(nPlots/2), 2, i);
        % Plot reference
        plot(tz, Nreal, 'k');
        hold all
        for j = 1:length(tz)
            Npt = zeros(1, nReps);
            for k = 1:nReps
                Npt = Nzhat{iPlots(i), k}(j); 
            end
            % Calculate mean and std deviation
            sig(i, j) = std(Npt);
            mu(i, j) = mean(Npt);
        end
        plot(tz, mu(i, :), 'r', tz, mu(i, :)-2*sig(i, :), 'g--', tz, mu(i, :)+2*sig(i, :), 'g--');
        xlabel('time');
        ylabel('population size');
        title(['Con-exp-con, ' xlab2{iPlots(i)}]);
        legend('N(t)', 'N(t) hat', 'bounds', 'location', 'best');
        grid;
        hold off
    end
    linkaxes(ax, 'xy');
            
    % Exponential function - bounds and estimate
    figure;
    Nreal1 = Nz1;
    % Get standard deviation bounds
    sig = zeros(nPlots, length(tz));
    mu = sig;
    for i = 1:nPlots
        ax(i) = subplot(ceil(nPlots/2), 2, i);
        % Plot reference
        plot(tz, Nreal1, 'k');
        hold all
        for j = 1:length(tz)
            Npt = zeros(1, nReps);
            for k = 1:nReps
                Npt = Nzhat1{iPlots(i), k}(j); 
            end
            % Calculate mean and std deviation
            sig(i, j) = std(Npt);
            mu(i, j) = mean(Npt);
        end
        plot(tz, mu(i, :), 'r', tz, mu(i, :)-2*sig(i, :), 'g--', tz, mu(i, :)+2*sig(i, :), 'g--');
        hold off
        xlabel('time');
        ylabel('population size');
        title(['Exponential, ' xlab2{iPlots(i)}]);
        legend('N(t)', 'N(t) hat', 'bounds', 'location', 'best');
        grid;
    end
    linkaxes(ax, 'xy');
    
    % Con-exp-con mean error at each time point for each scheme
    etraj = zeros(nRuns, length(tz));
    % Loop across no. schemes then time points then repetitions
    for i = 1:nRuns
        for j = 1:length(tz)
            ekk = zeros(1, nReps);
            for kk = 1:nReps
                % Error at jth time point for each repetition of ith scheme
                ekk(kk) = Nzhat{i, kk}(j) - Nreal(j);
            end
            % Mean traj error from all repetitions at each point
            etraj(i, j) = mean((ekk.^2));
        end
    end
    % Exponential mean error at each time point for each scheme
    etraj1 = zeros(nRuns, length(tz));
    % Loop across no. schemes then time points then repetitions
    for i = 1:nRuns
        for j = 1:length(tz)
            ekk = zeros(1, nReps);
            for kk = 1:nReps
                % Error at jth time point for each repetition of ith scheme
                ekk(kk) = Nzhat1{i, kk}(j) - Nreal1(j);
            end
            % Mean traj error from all repetitions at each point
            etraj1(i, j) = mean((ekk.^2));
        end
    end
    
    % Plot the mean square errors across the schemes
    figure;
    subplot(2, 1, 1);
    semilogy(tz, etraj);
    legend(xlab2, 'location', 'best');
    subplot(2, 1, 2);
    semilogy(tz, etraj1);
    legend(xlab2, 'location', 'best');
end


%% Plotting for the ranks and parameter estimates

% Plot ranks against schemes
figure;
subplot(2, 1, 1);
plot(dom, schemeMSE, 'o--', 'markersize', 8);
hold on
plot(dom, schemeMean, 's--', 'markersize', 8);
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
legend('MSE rank', 'Mean abs rank', 'location', 'best');
xlabel('no. sample times');
ylabel('rank - lowest is best');
title(['Rank of con-exp-con sampling schemes']);
subplot(2, 1, 2);
plot(dom, schemeMSE1, 'o--', 'markersize', 8);
hold on
plot(dom, schemeMean1, 's--', 'markersize', 8);
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
legend('MSE rank', 'Mean abs rank', 'location', 'best');
xlabel('no. sample times');
ylabel('rank - lowest is best');
title(['Rank of exponential sampling schemes']);

% Compare probability distributions from histograms
figure;
for i = 1:4
    subplot(2, 2, i);
    if i < 3
        h = histogram(err{i}(1, :), 20);
    else
        h = histogram(err1{i - 2}(1, :), 20);
    end
    h.Normalization = 'probability';
    hold all
    for j = 2:nRuns
        if i < 3
            h = histogram(err{i}(j, :), 20);
        else
            h = histogram(err1{i - 2}(j, :), 20);
        end
        h.Normalization = 'probability';
    end
    hold off
    grid;
    if i < 3
        xlabel(['x_' num2str(i) ' - x_' num2str(i) 'hat']);
        ylabel('probability');
        %legend('isochronous', 'heterochronous', 'location', 'best');
        legend(xlab2, 'location', 'best');
        title('Estimation error for con-exp-con');
    else
        xlabel(['x_' num2str(i - 2) ' - x_' num2str(i - 2) 'hat']);
        ylabel('probability');
        %legend('isochronous', 'heterochronous', 'location', 'best');
        legend(xlab2, 'location', 'best');
        title('Estimation error for exponential');
    end
end

% Plot the MSE for all parameters on a single figure by normalising to the
% maximum MSE scheme
figure;
subplot(2, 1, 1);
plot(dom, mseE{1}/max(mseE{1}), 'bo-', 'markersize', 8);
legcell = cell(1, numRV);
legcell{1} = 'x_1';
hold all
for i = 2:numRV
    plot(dom, mseE{i}/max(mseE{i}), 'o-', 'markersize', 8);
    legcell{i} = ['x_' num2str(i)];
end
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
xlabel('no. sample times');
legend(legcell, 'location', 'best');
title(['MSE ratio for con-exp-con at [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);
subplot(2, 1, 2);
plot(dom, mseE1{1}/max(mseE1{1}), 'bo-', 'markersize', 8);
legcell = cell(1, numRV);
legcell{1} = 'x_1';
hold all
for i = 2:numRV
    plot(dom, mseE1{i}/max(mseE1{i}), 'o-', 'markersize', 8);
    legcell{i} = ['x_' num2str(i)];
end
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
xlabel('no. sample times');
legend(legcell, 'location', 'best');
title(['MSE ratio for exponential at [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);


% Plot the absolute mean error for all parameters on a single figure by
% normalising to the maximum absolute error scheme
figure;
subplot(2, 1, 1);
plot(dom, abs(meanE{1})/max(abs(meanE{1})), 'bo-', 'markersize', 8);
legcell = cell(1, numRV);
legcell{1} = 'x_1';
hold all
for i = 2:numRV
    plot(dom, abs(meanE{i})/max(abs(meanE{i})), 'o-', 'markersize', 8);
    legcell{i} = ['x_' num2str(i)];
end
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
xlabel('no. sample times');
legend(legcell, 'location', 'best');
title(['Mean ratio for con-exp-con [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);
subplot(2, 1, 2);
plot(dom, abs(meanE1{1})/max(abs(meanE1{1})), 'bo-', 'markersize', 8);
legcell = cell(1, numRV);
legcell{1} = 'x_1';
hold all
for i = 2:numRV
    plot(dom, abs(meanE1{i})/max(abs(meanE1{i})), 'o-', 'markersize', 8);
    legcell{i} = ['x_' num2str(i)];
end
hold off
h = gca;
if length(h.XTickLabel) ~= length(xlab)
    h.XTick = 1:nRuns;
end
h.XTickLabel = xlab;
grid;
xlabel('no. sample times');
legend(legcell, 'location', 'best');
title(['Mean ratio for exponential [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);
