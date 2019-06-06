function [waveTimes, mea] = get_discharge_times(mea, varargin)
% Compute peak discharge times

p = inputParser;

addRequired(p, 'mea', @(x) isstruct(x) || strcmpi(class(x), 'matlab.io.MatFile'));
addParameter(p, 'verbose', true, @islogical);

parse(p, mea, varargin{:})
struct2var(p.Results)

if ~isstruct(mea)
	if ~exist(mea.Properties.Source, 'file')
		error('File not found')
	else
		mea = load(mea.Properties.Source); 
	end
end

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

padding = mea.Padding;
waveTimes = waveTimes - padding(1) * 1e3;  % Account for padding (in ms)

mea.waveTimes = waveTimes;
