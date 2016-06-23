% Cheat code to extract data from a figure 2paramExp1

% Assumes figure is open



D=get(gca,'Children'); %get the handle of the line object
XData=get(D,'XData');  %get the x data
YData=get(D,'YData');


figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    mxi = mx(i)*ones(1, M);
    %sxi = mx(i)*ones(1, M);
    histogram(mse_r(i, :), M/100);
    %plot(1:M, mse_r(i, :), 1:M, mxi, 'k', 1:M, mxi - 2*sxi, 'r', 1:M, mxi + 2*sxi, 'r')
    ylabel('freq')
    xlabel(['% (1 - x_' num2str(i) 'hat/x_' num2str(i) ')^2'])
    title(['J_r(x_' num2str(i) ') = ' num2str(mxi(1))]);
    %legend('relative mse', 'mean', 'mean - 2*std', 'mean+ 2*std', 'location', 'best');
    grid;
end

