function mea = mua_firing_rate(mea)
% Returns the firing rate as computed in Smith et al., 2016. Units are in
% spikes per second

events = zeros(size(mea.mua), 'single');
events(mea.event_inds) = true;
samplingRateMs = mea.SamplingRate / 1e3;  % samples per ms
windMs = 100;
window = windMs * samplingRateMs;  % number of samples to use in the window
mea.firingRate = ...
	smoothdata(events, 1, 'movmean', window) * ...
	mea.SamplingRate;  % convert to spikes per second 