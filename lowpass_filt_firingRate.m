function [frPct, Time] = lowpass_filt_firingRate(mea)

firingRate = mea.firingRate;
Time = mea.Time;
te = Time(end) - mea.Padding(1, 2); 
Time = Time(1:1e3:end);
firingRate = smoothdata(firingRate, 'movmedian', 1e3);  % smoothing over 30 ms window
firingRate = firingRate(1:1e3:end, :);  % Downsample to 30Hz

bpFilt = designfilt('lowpassfir','FilterOrder',150, ...
	'Passbandfrequency',.1, ...
	'Stopbandfrequency',1, ...
	'SampleRate',30);
frFilt = single(filtfilt(bpFilt, double(firingRate)));
% frFilt = smoothdata(frFilt(1:10:end, :));  % Downsample to 30 Hz
frPct = (frFilt - min(frFilt)) ./ (max(frFilt) - min(frFilt));
inds = logical((Time > 0) .* (Time < te));

frPct = frPct(inds, :); 
Time = Time(inds);

%%
figure(2), plot(Time, frPct + (1:size(frPct, 2)))