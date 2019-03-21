function [output, Time] = lowpass_filt_firingRate(mea, downsample)
% Tihs file is used to filter the firing rate of MEA data in preparations
% to make a video showing ictal wavefront spread. Data are downsampled to
% approximately 30 Hz and then low-pass filtered at 1 Hz. Finally, data are
% normalized to range [0, 1]. Only the ictal period of the data are
% returned.
% 
% Inputs
%    mea:			matfile containing properties firingRate, Time, Padding
%    downsample:	Downsample the data to this rate (Hz). 
%						Default: 30

if ~exist('downsample', 'var') || isempty(downsample)
	downsample = 30;
end

firingRate = mea.firingRate;
Time = mea.Time;
Time = Time();
te = Time(end) - mea.Padding(1, 2); 
SamplingRate = mea.SamplingRate;
skipfactor = round(SamplingRate / downsample);
Time = Time(1:skipfactor:end);
firingRate = firingRate(1:skipfactor:end, :);  % Downsample to 300Hz
SamplingRate = SamplingRate / skipfactor;
% firingRate = smoothdata(firingRate, 'movmean', SamplingRate);  % smoothing over 30 ms window

bpFilt = designfilt('lowpassfir', 'FilterOrder', 10, ...
	'CutoffFrequency', 1, ...
	'SampleRate', SamplingRate*1, ...
	'Window', 'hamming');
% frFilt = single(filtfilt(bpFilt, double(firingRate)));
frFilt = smoothdata(firingRate, 'gaussian', 5*SamplingRate);
% frFilt = smoothdata(frFilt(1:10:end, :));  % Downsample to 30 Hz
% frFilt = firingRate;
inds = logical((Time > 0) .* (Time < te));
output = frFilt(inds, :);

% frPct = (frFilt - min(frFilt)) ./ (max(frFilt) - min(frFilt));
% frPct = frPct(inds, :); 
Time = Time(inds);

%%
figure(2), plot(Time, output + .1*(1:size(output, 2)))
