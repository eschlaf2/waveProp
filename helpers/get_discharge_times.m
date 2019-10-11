function [waveTimes, mea] = get_discharge_times(mea, varargin)
% Compute peak discharge times

METHOD = mea.params.discharge_method;
MIN_FR = mea.params.min_fr;  % (method 1);
MIN_PD = mea.params.min_peak_distance;  % (method 1)
MIN_DEV = mea.params.min_lfp_deviance;  % [method 2; 2 (sd)]
MIN_ACTIVE = mea.params.min_active_electrodes;  % (method 2; default: 10)
WIN = mea.params.smoothing_win;  % [method 2; default: 10 (ms)];

if ~isstruct(mea)
	if ~exist(mea.Properties.Source, 'file')
		error('File not found')
	else
		mea = load(mea.Properties.Source); 
	end
end

switch METHOD
	case 1  % Based on peaks in firing rate

		try
			fr = mea.firingRate;
		catch ME
			if ~strcmp(ME.identifier, 'MATLAB:nonExistentField'), rethrow(ME); end
			disp('Computing firing rate.')
			[fr, mea] = mua_firing_rate(mea);
		end

		mask = mean(fr) >= MIN_FR;  % exclude electrodes that spike less than MIN_FR
		meanFr = mean(fr(:, mask), 2);

		[~, waveTimes] = findpeaks(meanFr, ...  % find peaks in mean firing rate
			mea.SamplingRate / 1e3, ...  % ... in ms 
			'minpeakprom', 100 * std(diff(meanFr)), ...  % ... use discrete peaks
			'minpeakdistance', MIN_PD);  % ... that are at least 100 ms apart
		waveTimes = waveTimes - mea.Padding(1) * 1e3;
		
	case 2  % based on peaks in LFP
		
		if ~isfield(mea, 'lfp'), [~, mea] = filter_mea(mea, 'lfp'); end
		[~, nCh] = size(mea.lfp);
		time = downsample(mea.Time(), mea.skipfactor) * 1e3;  % time in ms 

		mn = mean(mea.lfp(time < 0, :));  % preseizure mean
		sd = std(mea.lfp(time < 0, :));  % preseizure std
		normed_data = -(smoothdata((mea.lfp - mn) ./ sd));
		d_normed = [zeros(1, nCh); diff(normed_data)];
		peaks = ...  % Crossing 2sd below baseline:
			(d_normed > 0) & ...  % ... slope is positive on the left side
			(circshift(d_normed, -1) < 0) & ...  % ... and negative on the right side
			(normed_data > MIN_DEV);  % ... and is at least two sd from mean

		allpeaks = smoothdata(sum(peaks, 2), 'movmean', WIN) * WIN;  % show increases in activity within win
		waveTimes = time((allpeaks <= 10) & (circshift(allpeaks, -1) >= MIN_ACTIVE));  % Find where at least 10 electrodes peak within win of each other
end

mea.waveTimes = waveTimes;
