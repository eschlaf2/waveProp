function [event_inds, artefact_inds, mea, shape] = mua_events(mea)
% Takes a structure or matfile as input and computes event times. Events
% are defined as peaks in MUA more than EVENT_THRESH standard deviations from the
% baseline in the negative direction and at least MIN_DIST ms apart
% Uses global variable globs

EVENT_THRESH = mea.params.event_thresh;      % min peak distance from baseline (in negative direction) [sd]
MIN_DIST = mea.params.min_dist;          % min time between peaks [ms]
ARTEFACT_THRESH = mea.params.artefact_thresh;  % activity outside this range is considered artefact [sd] 

% Import fields
if (isfield(mea, 'mua') || isprop(mea, 'mua'))
	data = mea.mua;
else
	if ~isstruct(mea), mea = load(mea.Properties.Source); end
	mua = filter_mea(mea, 'mua');
	data = mua.mua;
	mea.mua = data;
end

time = mea.Time;
time = time();

% define artefacts as events that are more than ARTEFACT_THRESH sd from
% mean; remove them (set to nan)
artefacts = (abs(zscore(data)) > ARTEFACT_THRESH);  
data(artefacts) = nan;

intervalM = mea.SamplingRate / 1e3 * MIN_DIST;  % samples per MIN_DIST
if ~any(time < 0) 
	mn = 0; 
	sd = .01 * range(data);
else
	mn = mean(data(time < 0, :), 1, 'omitnan');     % get the mean of the preictal baseline
	sd = std(data(time < 0, :), 'omitnan');            % ... and the sd
end
data = (data - mn) ./ sd;                       % zscore based on preictal state

% find events on each channel
events = false(size(data));
for ch = 1:size(data, 2)
	temp = data(:, ch);
	if any(isnan(temp))
		i = find(isnan(temp), 1);
		while i
			i_next = find(~isnan(temp(i:end)), 1) + i;
			y = interp1([i - 1; i_next], [temp(i - 1); temp(i_next)], i:i_next - 1);
			temp(i:i_next - 1) = y;
			i = find(isnan(temp),1);
		end
	end
	if ~any(abs(zscore(temp)) > EVENT_THRESH)
		continue
	end
	[~, inds] = findpeaks(-(temp(~isnan(temp))), ...
		'minpeakdistance', intervalM, 'minpeakheight', EVENT_THRESH);  
	events(inds, ch) = true;
end

% store results to mea
shape = size(data);
event_inds = find(events);
artefact_inds = find(artefacts);

if nargout > 2
	mea.artefact_inds = artefact_inds;
	mea.event_inds = event_inds;
	mea.event_mat_size = shape;
end