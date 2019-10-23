% method = 'M';

% Params
rng default
pat = 'FHN'; seizure = 0; 
nTrials = 5;
sig_thresh = .05;
noise_levels = 2.^(5:-1:-2);  % SNR levels to test

% Setup
outfile = matfile(sprintf('testing%spinknoise%s%s_%d_%s', ...
	filesep, filesep, pat, seizure, method), 'Writable', true);
mea = load(sprintf('%s%s%s_Seizure%d_Neuroport_10_10.mat', ...
	pat, filesep, pat, seizure));
mea.params = init_params();
cn = dsp.ColoredNoise('InverseFrequencyPower', 2, ...
	'SamplesPerFrame', size(mea.Data, 1), ...
	'NumChannels', size(mea.Data, 2));
signal = mea.Data;
mask = zscore(signal) < 2;  % consider spikes as the signal
signal(mask) = nan;

rms =@(x) sqrt(mean(x.^2, 'omitnan'));
Asig = rms(signal);
noise =@(x, snr) x * (Asig / rms(x) / sqrt(snr));

switch upper(method(1))
	case 'D'
		if strcmpi(method, 'd') || strcmpi(method, 'd1')
			mea.params.T = 1;
			mea.params.delay_band = [1 50];
			method_full = 'delays_T01_fband1_50';
		else
			method_full = 'delays_T10_fband1_13';
		end
		method = 'delays';
	case 'M'
		method = 'maxdescent';
	case 'E'
		method = 'events';
end

% Get wave fits from full method
full = load(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pat, seizure), method_full);
full = full.(method_full);

% Use wavetimes from the full method on the subsets
mea.waveTimes = full.computeTimes;

% Convenience 
speed =@(v) sqrt(sum(v.^2));
nWaves = numel(mea.waveTimes);

% Reset mea to original after each trial
original_mea = mea;

for snr = noise_levels
	
	% Initialize stats
	stats.detections = zeros(nTrials, 1);
	stats.detection_rate = nan(nTrials, 1);
	[stats.dtheta, stats.dspeed] = deal(nan(nTrials, nWaves));
	
	for trial = 1:nTrials  % For each trial
		
		% Add noise
		mea.Data = mea.Data + noise(cn(), snr);
		
		% Compute the wave fits and save to outfile
		fits = wave_prop(mea, method);
		fits = rmfield(fits, {'data', 'position'});
		fieldname = checkname(sprintf('snr%02g_%d', snr, trial));
		outfile.(fieldname) = fits;  % save details
		
		% Store stats of fits compared to full
		mask = fits.p < sig_thresh & full.p < sig_thresh;
		stats.detections(trial) = sum(mask);
		stats.detection_rate(trial) = ...
			sum(mask) ./ sum(full.p < sig_thresh);
		stats.dtheta(trial, mask) = exp(1j*fits.Z(mask)) - exp(1j*full.Z(mask));
		stats.dspeed(trial, mask) = speed(fits.V(:, mask)) - speed(full.V(:, mask));
		mea = original_mea;
	end
	
	% Save stats of fits to outfile
	outfile.(checkname(sprintf('snr%02d_stats', snr))) = stats;
end