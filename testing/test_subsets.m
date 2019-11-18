% fname = 'CUCX4/CUCX4_Seizure2_Neuroport_10_10.mat'; method = 'M';

% Params
rng default
nTrials = 100;
sig_thresh = .05;
subsets = [4 9 16 25 36 49 64];

% Setup
mea = load(fname);
mea.params = init_mea_params();
[~, fname, ~] = fileparts(fname);
outfile = matfile(sprintf('testing%ssubsets%s%s_%s', ...
	filesep, filesep, fname, method), 'Writable', true);

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

if ~exist('method_full', 'var'); method_full = method; end

% Get wave fits from full method
full = load(sprintf('%s_wave_prop.mat', fname), method_full);
full = full.(method_full);

% Use wavetimes from the full method on the subsets
mea.waveTimes = full.computeTimes;

% Convenience 
speed =@(v) sqrt(sum(v.^2));
nWaves = numel(mea.waveTimes);

% Reset mea to original after each trial
original_mea = mea;

for sub = subsets
	
	% Initialize stats
	stats.detections = zeros(nTrials, 1);
	stats.detection_rate = nan(nTrials, 1);
	[stats.dtheta, stats.dspeed] = deal(nan(nTrials, nWaves));
	
	for trial = 1:nTrials  % For each trial
		
		% select a subset of channels
		issub = keep_random_channels(mea, sub);
		if ~issub, continue, end  % exit if the subset is the size of the number of channels
		
		% Compute the wave fits and save to outfile
		fits = wave_prop(mea, method);
		fits = rmfield(fits, {'data', 'position'});
		fieldname = sprintf('sub%02d_%d', sub, trial);
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
	outfile.(sprintf('sub%02d_stats', sub)) = stats;
end