% Function to resample uniformly the last L values of data
function [tx, Nhatx, lamhatx] = makeUniformSamples(tn, L, F, Nhat, lamhat)

% Assumptions and modifications
% - only last L values are important

% Check lengths
lentn = length(tn);
if L > lentn
    error('Specified L is too large');
end

% Get resample time for last L values of tn
idchange = (lentn - (L + 1)):lentn;
tnchange = tn(idchange);
dt = mean(diff(tnchange))/F;

% Construct new time vector and the tn segment for interpolating
tx = tnchange(end):-dt:tnchange(1);
tx = tx(end:-1:1);
lentx = length(tx);
if lentx > 10^6
    warning('Mat:length', 'Length of resampled vector is too long');
    asignin('base', 'lentx', lentx);
end
tx = tx';

% Linearly interpolate for parameters of interest
Nhatx = nakeinterp1(tnchange, Nhat(idchange), tx);
lamhatx = nakeinterp1(tnchange, lamhat(idchange), tx);