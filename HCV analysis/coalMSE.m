% Function to calculate MSE from coalescent data
function [statSet, em, tm, setData] = coalMSE(tn, N, Nhat)

% Assumptions and Modifications
% - tm is not uniformly sampled, it is just a non-duplicate form of tn
% - added running MSE to see how it progresses across samples
% - removed sets of size zero or 1 which lead to NaNs in calc %<----------

% Resample the time series by first removing any initial times and deciding
% a sample interval, delt and get raw error, en
tn = tn - tn(1);
%delt = mean(diff(tn))/100; % may have issues in coalescent due to diff t scales
dtn = diff(tn);
dtn = dtn(dtn > 0);
delt = min(dtn)/10;
nSamps = range(tn)/delt;
en = N - Nhat;

% Remove duplicate points - the +1 in the id chooses the latter of point
dtn = diff(tn);
id = find(dtn == 0) + 1; 
lentn = length(tn);

% Obtain altered data vectors with duplicate time indices removed
idset = 1:lentn;
trunc = setdiff(idset', id);
tm = tn(trunc);
em = en(trunc);

% Separate samples in sets of Tsetsize
Tsetsize = 20000;
nSetSamps = floor(nSamps/Tsetsize);
tSetSamps = range(tn)/nSetSamps;

% Linear interpolation variable declaration for en, tlim(1) set for no bias
sum_em = 0;
sum_emSq = 0;
sum_len = 0;
tlim = -ones(1, nSetSamps+1);
tlim(1) = 0;

% Store square error across samples and no. samples in set
em_run_mse = zeros(1, nSetSamps);
run_samplen = zeros(1, nSetSamps);

% Loop across sets iteratively calculating statistics
for i = 2:nSetSamps+1
    % Define time limits and obtain relevant section of data
    tlim(i) = tSetSamps*(i-1);
    idtemp = find(tm >= tlim(i-1) & tm < tlim(i));
    
    % Account for empty ttemp if the sampling is small
    if ~isempty(idtemp) && length(idtemp) > 1
        ttemp = tm(idtemp);
        etemp = em(idtemp);
        
        % Obtain sample times and interpolate em data to the sampled times
        tsamp = ttemp(1):delt:ttemp(end);
        esamp = nakeinterp1(ttemp, etemp, tsamp');
        
        lenSamp = length(esamp);
        if lenSamp ~= length(tsamp)
            assignin('base', 'tsamp', tsamp);
            assignin('base', 'eSamp', esamp);
            error(['The sampling of the error curve failed at i = ' num2str(i)]);
        end
        
        % Obtain iterative sums for statistics
        esampSq = esamp.*esamp;
        sum_em = sum_em + sum(esamp);
        sum_emSq = sum_emSq + sum(esampSq);
        sum_len = sum_len + lenSamp;
        
        % Obtain running MSE
        em_run_mse(i) = sum_emSq/sum_len;
        run_samplen(i) = lenSamp;
    end
end

% Check that the correct number of points were obtained and calculate means
if sum_len > nSamps 
    assignin('base', 'tlim', tlim);
    error(['Inconsistent sample size: [nSamps sum_len] = ' [num2str(nSamps)...
        ' ' num2str(sum_len)]]);
else
    em_mean = sum_em/sum_len;
    em_mse = sum_emSq/sum_len;
    em_var = em_mse - em_mean^2;
end

% Assign stats set with indicator string to indicate data order
statSet.order = {'mean', 'var', 'mse'};
statSet.val = [em_mean em_var em_mse];
setData{1} = em_run_mse; 
setData{2} = run_samplen;