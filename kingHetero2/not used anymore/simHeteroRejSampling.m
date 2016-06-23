% Function to simulate a heterochronous coalescent using rejection sampling
function events = simHeteroRejSampling(fn)

% Assumptions and modifications
% - also calculate no. lineages since it isn't a simple death process
% - svec and nvec are sampling times and no. lineages introduced
% - Lmax assumed to be in ascending order so it is indexed by nLinCurr

% Total number of lineages and samples, extract from input
nData = fn.nData;
n = nData + 1;
nvec = fn.nvec;
svec = fn.svec;
nSampTimes = length(svec);
Lmax = fn.Lmax;

% Ensure distince sample times and the first sample time is 0
if length(unique(svec)) ~= nSampTimes
    error('The sample times are not distinct');
else
    if svec(1) ~= 0
        error('First sample time should match an isochronous case');
    end
end

% Variables for coalescent times and total no. lineages at all event times
tcoal = zeros(1, n);
nLin = zeros(1, nData + nSampTimes);
tLin = zeros(1, nData + nSampTimes);

% Logging indices for Gillespie type algorithm with deterministic events
I = 1; % for coalescent events but n-1 coalescents so tcoal(1) = svec(1)
J = 2; % for sampling events, as J = 1 is start point at t = 0
t = 0;
nLin(1) = nvec(1);
nLinCurr = nLin(1);
ind = 2; % counts for all events, first event is start
sampleTrue = 0;

% Simulate heterochronous coalescent process with thinning and treatment of
% discontinuous, deterministic sample introduction
while(I <= nData)
    % Ensure no. lineages never goes to zero before all samples coalesce
    if (nLinCurr == 1 || sampleTrue) && J <= nSampTimes
        % Injection of new lineages
        nLinCurr = nLinCurr + nvec(J);
        nLin(ind) = nLinCurr;
        
        % Sampling event - discontinuously update time and tcurr = t
        tLin(ind) = svec(J);
        t = svec(J);
        sampleTrue = 0;
        J = J + 1;
        ind = ind + 1;      
    else
        % Generate a Poisson homogeneous interval from last time
        U = rand;
        t = t -log(U)/Lmax(nLinCurr);
        
        % Account for J value after last J+1 update to precent vec overflow
        if J > nSampTimes
            nextSampTime = inf;
        else
            nextSampTime = svec(J);
        end
        
        % Decide if event is coalescent or sampling
        if t <= nextSampTime
            % Calculate rate at current time
            bfac = 0.5*nLinCurr*(nLinCurr-1);
            [~, lamt] = getTimeVaryingN(fn, t, bfac, fn.param');
            
            % Rejection sample by calculating rate
            U = rand;
            if U <= lamt/Lmax(nLinCurr)
                % An event has occurred so take data, ensure tcoal(1) = 0
                I = I + 1;
                tcoal(I) = t;
                
                % Lineages fall by 1 due to coalescent
                nLin(ind) = nLinCurr -1;
                nLinCurr = nLinCurr - 1;
                tLin(ind) = t;
                ind = ind + 1;
            end
            
        else
            % Switch to next sample
            sampleTrue = 1;
        end
    end
    
%     disp(num2str(I));
%     disp(num2str(ind));
%     disp(num2str(nLinCurr));
%     disp(num2str(lamt));
end

% Check for consistency in event times
tcheck = sort([tcoal svec(2:end)]);
if ~all(tcheck == tLin)
    events.tcheck = tcheck;
    disp('The times are inconsistent');
end

% A total event time vector containing coalescents and samples - take
% svec(2:end) as svec(1) = tcoal(1) = 0
events.tLin = tLin;
events.nLin = nLin;
events.tcoal = tcoal;