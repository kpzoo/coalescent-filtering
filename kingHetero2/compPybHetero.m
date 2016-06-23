% Script to compare batch runs on pybus N(t) that are
% heterochronously sampled: con-exp-con parametrised as in HCV

% Assumptions and Modifications
% - 

clearvars
clc
close all

%% Set loop parameters and run simulation and filtering for exponential

% Initialise total number of samples and use factors to get nSampTimes
n = 100; 
zz = 1:n;
nSampTimes = zz(~(rem(n, zz)));
nSampTimes = nSampTimes(1:end-1); % don't want nSampTimes = n case
sampSize = n./nSampTimes;

% Define no. runs and repetitions and maxSamp time and filter dimension
nReps = 2000;
nRuns = length(nSampTimes);
maxSamp = 100;
numRV = 4;
mi = 15*ones(1, numRV);

% Variables to store outputs - each row is a different sampling scheme
ss = zeros(nRuns, nReps);
neffM = ss; cap = ss; lamM = ss; lamV = ss;
x1 = ss; x1hat = ss; x2 = ss; x2hat = ss; x3 = ss; x3hat = ss; x4 = ss; x4hat = ss;

% Main loop to repeat many iterations of the same sampling and then also
% loop across different sampling schemes
for ii = 1:nRuns
    for jj = 1:nReps
        % Main simulation of coalescent with filtering
        [neffM(ii, jj), cap(ii, jj), lamM(ii, jj), lamV(ii, jj), est] = pybHeteroFn(n, nSampTimes(ii), maxSamp, mi);
        % Parameters of interest
        x1(ii, jj) = est.x(1); x2(ii, jj) = est.x(2);
        x3(ii, jj) = est.x(3); x4(ii, jj) = est.x(4);
        x1hat(ii, jj) = est.xhat(1); x2hat(ii, jj) = est.xhat(2);
        x3hat(ii, jj) = est.xhat(3); x4hat(ii, jj) = est.xhat(4);
    end
    disp(['Finished sampling scheme ' num2str(ii) ' of ' num2str(nRuns)]);
end

% Save the data
save(['pybSet' num2str(n)]);

% Estimation statistics
e1 = (x1 - x1hat)./x1; e2 = (x2 - x2hat)./x2;
e3 = (x3 - x3hat)./x3; e4 = (x4 - x4hat)./x4;
e1sq = e1.^2; e2sq = e2.^2; e3sq = e3.^2; e4sq = e4.^2;

% Collect data into more manageable and extendable format
err{1} = e1; err{2} = e2; err{3} = e3; err{4} = e4;
meanE = cell(1, numRV);
mseE = meanE;
for i = 1:numRV
    esq = err{i}.^2;
    meanE{i} = mean(err{i}, 2);
    mseE{i} = mean(esq, 2);
end

% Mean and square error values for each scheme
neff = mean(neffM, 2);
C = mean(cap, 2);
lamMean = mean(lamM, 2); lamVar = mean(lamV, 2);

% Performs simple ranking procedure - assign value nReps to worst MSE or
% mean error and 1 to lowest
MSErnks = zeros(numRV, nRuns);
Meanrnks = MSErnks;
for i = 1:numRV
    MSErnks(i, :) = tiedrank(mseE{i});
    Meanrnks(i, :) = tiedrank(abs(meanE{i}));
end

% Then get overall rank based on summing ranks of each scheme for each xi
schemeMSE = tiedrank(sum(MSErnks));
schemeMean = tiedrank(sum(Meanrnks));

%% Plotting and data visualisation

% Control the x axis labelling and form
dom = 1:nRuns;
xlab = cell(1, nRuns);
for i = dom
    xlab{i} = num2str(nSampTimes(i));
end

% Plot ranks against schemes
figure;
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
title(['Rank of various sampling schemes']);

% Compare histograms for each variable and sample scheme
freq = cell(numRV, nRuns);
val = freq;
figure;
for i = 1:numRV
    subplot(numRV/2, numRV/2, i);
    [freq{i, 1}, val{i, 1}] = hist(err{i}(1, :), 20);
    plot(val{i, 1}, freq{i, 1}, 'bo-');
    hold all
    for j = 2:nRuns
        [freq{i, j}, val{i, j}] = hist(err{i}(j, :), 20);
        plot(val{i, j}, freq{i, j});
    end
    hold off
    grid;
    xlabel(['x_' num2str(i) ' - x_' num2str(i) 'hat']);
    ylabel('freq');
    legend('isochronous', 'heterochronous', 'location', 'best');
    title('Estimation error histograms');
end

% Compare probability distributions from histograms
figure;
for i = 1:numRV
    subplot(numRV/2, numRV/2, i);
    h = histogram(err{i}(1, :), 20);
    h.Normalization = 'probability';
    hold all
    for j = 2:nRuns
        h = histogram(err{i}(j, :), 20);
        h.Normalization = 'probability';
    end
    hold off
    grid;
    xlabel(['x_' num2str(i) ' - x_' num2str(i) 'hat']);
    ylabel('probability');
    legend('isochronous', 'heterochronous', 'location', 'best');
    title('Estimation error histograms');
end

% Plot the MSE for all parameters on a single figure by normalising to the
% maximum MSE scheme
figure;
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
title(['MSE ratio at [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);

figure;
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
title(['Mean ratio at [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);

