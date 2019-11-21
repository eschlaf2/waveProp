function [fr, mea] = mua_firing_rate(mea)
% Returns the firing rate as computed in Smith et al., 2016. Units are in
% spikes per second

% Load necessary variables (faster with matfile objects)
if ~isstruct(mea), mea = load(mea.Properties.Source); end
if ~isfield(mea, 'event_inds'), [~, ~, mea, ~] = mua_events(mea); end

event_inds = mea.event_inds;
samplingRate = mea.SamplingRate;
Time = mea.Time;
T = numel(Time());
NCh = size(mea.Data, 2);
if any(strcmpi(fieldnames(mea), 'badchannels'))
	NCh = NCh - numel(mea.BadChannels);
end
dims = [T, NCh];

% Convert event_inds into binned spike increments
events = zeros(dims, 'single');
events(event_inds) = 1;
samplingRateMs = samplingRate / 1e3;  % samples per ms
window = mea.params.fr_window * samplingRateMs;  % number of samples to use in the window

% Compute spike rate
fr = smoothdata(events, 1, 'movmean', window) * ...
	samplingRate;  % convert to spikes per second 

mea.firingRate = fr;
