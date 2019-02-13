function [mea] = mua_events(mea)
% Takes a structure or matfile as input and computes event times. Events
% are defined as peaks in MUA more than EVENT_THRESH standard deviations from the
% baseline in the negative direction and at least MIN_DIST ms apart

EVENT_THRESH = 4;      % min peak distance from baseline (in negative direction) [sd]
MIN_DIST = 1;          % min time between peaks [ms]
ARTEFACT_THRESH = 16;  % activity outside this range is considered artefact [sd] 

% Import fields
data = mea.mua;
time = mea.Time;

intervalM = mea.SamplingRate / 1e3 * MIN_DIST;  % samples per MIN_DIST
mn = mean(data(time < 0, :), 1);   % get the mean of the baseline
sd = std(data);                              % and the overall sd
temp = (data - mn) ./ sd;                    % zscore based on baseline mean

% define artefacts as events that are more than ARTEFACT_THRESH sd from
% baseline; remove them (set to nan)
artefacts = (abs(temp) > ARTEFACT_THRESH);  
data(artefacts) = nan;

% find events on each channel
events = false(size(temp));
for ch = 1:size(temp, 2)
	temp = data(:, ch);
	[~, inds] = findpeaks(-zscore(temp(~isnan(temp))), ...
		'minpeakdistance', intervalM, 'minpeakheight', EVENT_THRESH);  
	events(inds, ch) = true;
end

% store results to mea
mea.artefact_inds = find(artefacts);
mea.event_inds = find(events);