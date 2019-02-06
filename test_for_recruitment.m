function mea = test_for_recruitment(mea, PLOT)
% Test to see if the electrode array was recruited to the seizure as
% described in ?Schevon, Catherine A., et al. ?Evidence of an Inhibitory
% Restraint of Seizure Activity in Humans.? Nature Communications, vol. 3,
% no. 1, Nature Publishing Group, Jan. 2012, p. 1060,
% doi:10.1038/ncomms2056.
% 
% Uses weighted phase locking index

if ~exist('PLOT', 'var')
	PLOT = true;
end

%% Downsample lfp
ds_freq = 1e3;  % Hz
ds_step = floor(mea.SamplingRate / ds_freq);
ds_freq = mea.SamplingRate / ds_step;
lfp = mea.lfp(1:ds_step:end, :);
T = mea.Time(1:ds_step:end);

%% Set parameters for coherence calculations
TW = 20;
ntapers = 2*TW-1;
params.Fs = ds_freq;
params.tapers = [TW, ntapers];
params.pad = 0;
params.fpass = [2 30];
params.err = [1 .005];

%% Compute the multitaper complex cross spectrum at 160 ms intervals
window = 1;  % length of window (s)
winSamples = window * ds_freq;  % samples per window
winstep = floor(160e-3 * ds_freq);  % step through the data in 160 ms increments
nSamples = size(lfp, 1) - winSamples + 1;  % number of samples for which to compute coherence

% initialize vector to store coherence
wpli_mean = zeros(size(1:winstep:nSamples));  

i = 0;  % increment by one
for ti = 1:winstep:nSamples  % step through time at winstep increments
	i = i + 1;
	
	% Compute the complex cross spectrum
	[Sc,~,~,~,~,f] = CrossSpecMatc(...
		lfp(ti : ti + winSamples-1, :), window, params);
	
	% Compute the weighted phase locking index for each channel pair
	wpli = squeeze(abs(mean(imag(Sc), 1)) ./ mean(abs(imag(Sc)), 1));
	
	% take the mean of the squared wpli values
	wpli_mean(i) = mean(wpli(:).^2, 'omitnan');
end

if PLOT
	figure(4); fullwidth()
	plot(T(1:winstep:ti), wpli_mean);
	title([strrep(mea.Name, '_', ' '), ' Coherence'])
	ylabel('Mean square WPLI')
	xlabel('Time (s)')
	axis tight
end
	
mea.wpli = wpli_mean;
mea.wpliF = f;

