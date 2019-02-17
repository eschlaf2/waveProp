function mea = mua_firing_rate(mea)
% Returns the firing rate as computed in Smith et al., 2016. Units are in
% spikes per second

% Load necessary variables (faster with matfile objects)
event_inds = mea.event_inds;
dims = size(mea.mua);
samplingRate = mea.SamplingRate;

% Convert event_inds into binned spike increments
events = zeros(dims, 'single');
events(event_inds) = 1;
samplingRateMs = samplingRate / 1e3;  % samples per ms
windMs = 100;
window = windMs * samplingRateMs;  % number of samples to use in the window

% Compute spike rate
fr = smoothdata(events, 1, 'movmean', window) * ...
	samplingRate;  % convert to spikes per second 

mea.firingRate = fr;