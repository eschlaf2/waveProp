function [waveTimes, mea] = get_discharge_times(mea, varargin)
% Compute peak discharge times

p = inputParser;

addRequired(p, 'mea', @(x) isstruct(x) || strcmpi(class(x), 'matlab.io.MatFile'));
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'method', 2, @(x) any(x == [1 2]));

parse(p, mea, varargin{:})
struct2var(p.Results)

if ~isstruct(mea)
	if ~exist(mea.Properties.Source, 'file')
		error('File not found')
	else
		mea = load(mea.Properties.Source); 
	end
end

switch method
	case 1  % Based on peaks in firing rate

		try
			fr = mea.firingRate;
		catch ME
			if ~strcmp(ME.identifier, 'MATLAB:nonExistentField'), rethrow(ME); end
			disp('Computing firing rate.')
			[fr, mea] = mua_firing_rate(mea);
		end

		mask = mean(fr) >= 1/60;  % exclude channels with mean firing rate less than one spike per minute (Liou et al., 2018) ?\cite{Liou2018a}
		meanFr = mean(fr(:, mask), 2);

		[~, waveTimes] = findpeaks(meanFr, ...  % find peaks in mean firing rate
			mea.SamplingRate / 1e3, ...  % ... in ms 
			'minpeakprom', 100 * std(diff(meanFr)), ...  % ... use discrete peaks
			'minpeakdistance', 100);  % ... that are at least 100 ms apart
		time = mea.Time() * 1e3;
		waveTimes = time(waveTimes);
		
	case 2  % based on peaks in LFP
		
		if ~isfield(mea, 'lfp'), [~, mea] = filter_mea(mea, 'lfp'); end
		[~, nCh] = size(mea.lfp);
		time = downsample(mea.Time(), mea.skipfactor) * 1e3;  % time in ms 

		mn = mean(mea.lfp(time < 0, :)); 
		sd = std(mea.lfp(time < 0, :));
		normed_data = -(smoothdata((mea.lfp - mn) ./ sd));
		d_normed = [zeros(1, nCh); diff(normed_data)];
		peaks = ...  % Crossing 2sd below baseline:
			(d_normed > 0) & ...  % ... slope is positive on the left side
			(circshift(d_normed, -1) < 0) & ...  % ... and negative on the right side
			(normed_data > 2);  % ... and is at least two sd from mean


		win = 10;  % ms windows (lfp is already downsampled to 1 kHz)
		peaksF = diff(smoothdata(peaks, 'movmean', win));
		peaksF = (peaksF > 0) & (circshift(peaksF, -1) <= 0);
		allpeaks = smoothdata(sum(peaksF, 2), 'movmean', win) * win;  % show increases in activity within win
		waveTimes = time((allpeaks <= 10) & (circshift(allpeaks, -1) > 10));  % Find where at least 10 electrodes peak within win of each other
end

mea.waveTimes = waveTimes;
