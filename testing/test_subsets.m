% pat = 'CUCX4'; seizure = 2; method = 'M';

rng default
nTrials = 5;
sig_thresh = .05;
outfile = matfile(sprintf('testing%ssubsets%s%s_%d_%s', ...
	filesep, filesep, pat, seizure, method), 'Writable', true);

% Setup
full = load(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pat, seizure), method);
full = full.(method);
mea = load(sprintf('%s%s%s_Seizure%d_Neuroport_10_10.mat', ...
	pat, filesep, pat, seizure));

if strcmpi(method, 'D')
	mea.params = init_params('T', 1, 'delay_band', [1 50]); 
else
	mea.params = init_params();
end

switch upper(method(1))
	case 'D'
		method = 'delays';
	case 'M'
		method = 'maxdescent';
	case 'E'
		method = 'events';
end

mea.waveTimes = full.computeTimes;
nWaves = numel(mea.waveTimes);

% Convenience functions
speed =@(v) sqrt(sum(v.^2));

% Subsets of electrodes
subsets = [4 9 16 25 36 49 64];
original_mea = mea;

for sub = subsets
	% Initialize stats
	stats.detections = zeros(nTrials, 1);
	stats.detection_rate = nan(nTrials, 1);
	[stats.dtheta, stats.dspeed] = deal(nan(nTrials, nWaves));
	for trial = 1:nTrials
		issub = keep_random_channels(mea, sub);
		if ~issub, continue, end
		fits = wave_prop(mea, method);
		fits = rmfield(fits, {'data', 'position'});
		fieldname = sprintf('sub%02d_%d', sub, trial);
		outfile.(fieldname) = fits;  % save details
		
		% Save stats
		mask = fits.p < sig_thresh;
		stats.detections(trial) = sum(mask);
		stats.detection_rate(trial) = ...
			sum(mask) ./ sum(full.p < sig_thresh);
		stats.dtheta(trial, mask) = exp(1j*fits.Z(mask)) - exp(1j*full.Z(mask));
		stats.dspeed(trial, mask) = speed(fits.V(:, mask)) - speed(full.V(:, mask));
		mea = original_mea;
	end
	outfile.(sprintf('sub%02d_stats', sub)) = stats;
end