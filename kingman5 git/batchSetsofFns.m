% Script to batch several N(t) functions
clear all
clc 
close all

% Set key parameter of function ids, mi and no. repetitions M
fnIDs = [4 5 2];
lenfn = length(fnIDs);
M = 1000;
mi = 10*ones(1, 4);
Tavg = zeros(1, lenfn);
numRV = length(mi);

% Loop across functions
for iii = 1:lenfn
    [mse_t, param, paramhat, n, M, tIter, fn, C, Tavg(iii)] = batchMV2Fn(fnIDs(iii), M, mi);
    disp('***************************************************************');
    disp(['Finished function ' num2str(iii) ' of ' num2str(M)]);
    disp(['Avg time for ' num2str(M) ' iterations is ' num2str(Tavg(iii))]);
    disp('***************************************************************');
    % Save important bits of data
    save(['batch_' num2str(numRV) '_' num2str(fnIDs(iii)) '_' num2str(M)], 'mse_t', 'param', 'paramhat',...
        'n', 'M', 'tIter', 'fn', 'C', 'numRV');
end