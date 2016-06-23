% Script to compare batch runs on exponential N(t) that are
% heterochronously sampled

% Assumptions and Modifications
% - 

clearvars
clc
close all

%% Set loop parameters and run simulation and filtering for exponential

% Initialise total number of samples and use factors to get nSampTimes
n = 50; 
zz = 1:n;
nSampTimes = zz(~(rem(n, zz)));
nSampTimes = nSampTimes(1:end-1); % don't want nSampTimes = n case
sampSize = n./nSampTimes;

% Define no. runs and repetitions and maxSamp time and filter dimension
nReps = 5000;
nRuns = length(nSampTimes);
maxSamp = 10;
mi = [30 30];

% Variables to store outputs - each row is a different sampling scheme
ss = zeros(nRuns, nReps);
neffM = ss;
cap = ss;
lamM = ss;
lamV = ss;
x1 = ss;
x1hat = ss;
x2 = ss;
x2hat = ss;

% Main loop to repeat many iterations of the same sampling and then also
% loop across different sampling schemes
for ii = 1:nRuns
    for jj = 1:nReps
        % Main simulation of coalescent with filtering
        [neffM(ii, jj), cap(ii, jj), lamM(ii, jj), lamV(ii, jj), est] = expHeteroFn(n, nSampTimes(ii), maxSamp, mi);
        % Parameters of interest
        x1(ii, jj) = est.x(1);
        x2(ii, jj) = est.x(2);
        x1hat(ii, jj) = est.xhat(1);
        x2hat(ii, jj) = est.xhat(2);
    end
    disp(['Finished sampling scheme ' num2str(ii) ' of ' num2str(nRuns)]);
end

% Save the data
save(['expSet' num2str(n)]);

% Estimation statistics
e1 = (x1 - x1hat)./x1;
e2 = (x2 - x2hat)./x2;
e1sq = e1.^2;
e2sq = e2.^2;

% Mean and square error values for each scheme
meanE1 = mean(e1, 2);
meanE2 = mean(e2, 2);
mseE1 = mean(e1sq, 2);
mseE2 = mean(e2sq, 2);
neff = mean(neffM, 2);
C = mean(cap, 2);
lamMean = mean(lamM, 2);
lamVar = mean(lamV, 2);


%% Plotting and data visualisation

% Control the x axis labelling and form
dom = 1:nRuns;
xlab = cell(1, nRuns);
for i = dom
    xlab{i} = num2str(nSampTimes(i));
end

% Compare histograms for each variable and sample scheme
freq1 = cell(1, nRuns);
val1 = freq1;
val2 = val1;
freq2 = freq1;

figure;
subplot(1, 2, 1)
[freq1{1}, val1{1}] = hist(e1(1, :), 20);
plot(val1{1}, freq1{1}, 'bo-');
hold all
for i = 2:nRuns
    [freq1{i}, val1{i}] = hist(e1(i, :), 20);
    plot(val1{i}, freq1{i});
end
hold off 
grid;
xlabel('x_1 - x_1 hat');
ylabel('freq');
legend('isochronous', 'heterochronous', 'location', 'best');
title('Estimation error histograms');
subplot(1, 2, 2)
[freq2{1}, val2{1}] = hist(e2(1, :), 20);
plot(val2{1}, freq2{1}, 'bo-');
hold all
for i = 2:nRuns
    [freq2{i}, val2{i}] = hist(e2(i, :), 20);
    plot(val2{i}, freq2{i});
end
hold off 
grid;
xlabel('x_2 - x_2 hat');
ylabel('freq');
legend('isochronous', 'heterochronous', 'location', 'best');
title('Estimation error histograms');

% Compare probability distributions from histograms
figure;
subplot(1, 2, 1);
h = histogram(e1(1, :), 20);
h.Normalization = 'probability';
hold all
for i = 2:nRuns
    h = histogram(e1(i, :), 20);
    h.Normalization = 'probability';
end
hold off 
grid;
xlabel('x_1 - x_1 hat');
ylabel('probability');
legend('isochronous', 'heterochronous', 'location', 'best');
title('Estimation error histograms');
subplot(1, 2, 2);
h = histogram(e2(1, :), 20);
h.Normalization = 'probability';
hold all
for i = 2:nRuns
    h = histogram(e2(i, :), 20);
    h.Normalization = 'probability';
end
hold off 
grid;
xlabel('x_2 - x_2 hat');
ylabel('probability');
legend('isochronous', 'heterochronous', 'location', 'best');
title('Estimation error histograms');

% Plot the MSE for both parameters on a single figure
figure;
[hAx, hLine1, hLine2] = plotyy(dom, 100*mseE1, dom, 100*mseE2);
ylabel(hAx(1), 'perc mse x_1') % left y-axis
ylabel(hAx(2), 'perc mse x_2') % right y-axis
hLine1.LineStyle = '--';
hLine2.LineStyle = '--';
hLine1.Marker = 'o';
hLine2.Marker = 'o';
hLine1.MarkerSize = 8;
hLine2.MarkerSize = 8;
%hLine1.LineWidth = 2;
%hLine2.LineWidth = 2;
if length(hAx(1).XTickLabel) ~= length(xlab)
    hAx(1).XTick = 1:nRuns;
    hAx(2).XTick = 1:nRuns;
end
hAx(1).XTickLabel = xlab;
hAx(2).XTickLabel = xlab;
grid;
xlabel('no. sample times');
legend('%mse(x_1)', '%mse(x_2)', 'location', 'best');
title(['Relative MSE at [n mi nreps] = [' [num2str(n) ' ' num2str(mi(1)) ' ' num2str(nReps)] ']']);

