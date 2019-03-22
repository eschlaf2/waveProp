function waveTimes = get_discharge_times(mea)

fr = mea.firingRate;
mask = mean(fr) >= 1/60;  % exclude channels with mean firing rate less than one spike per minute (Liou et al., 2018) ?\cite{Liou2018a}
meanFr = mean(fr(:, mask), 2);

[~, waveTimes] = findpeaks(meanFr, ...  % find peaks in mean firing rate
	mea.SamplingRate / 1e3, ...  % ... in ms 
	'minpeakprom', 100 * std(diff(meanFr)), ...  % ... use discrete peaks
	'minpeakdistance', 100);  % ... peaks should be at least 100 ms apart

padding = mea.Padding;
waveTimes = waveTimes - padding(1) * 1e3;  % Account for padding (in ms)

mea.waveTimes = waveTimes;
