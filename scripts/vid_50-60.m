% Creates a video of the electrode array for the matfile variable mea,
% which should be in the workspace already.
% Example:
%     mea = matfile('c5/c5_Seizure1_Neuroport.mat');
%     electrode_vid_full;

mea = load('c5/c5_Seizure1_Neuroport_10_10.mat');  % load data
mea = exclude_channels(mea);  % remove non-participating channels
P = mea.Position;  % get electrode coordinates
P(mea.BadChannels, :) = [];  % ... and remove bad channels
Time = mea.Time;  % Get sampling times
Time = Time();  % ... (and convert from function to array if necessary)
te = Time(end) - mea.Padding(1, 2);  % Find when seizure ends (so that video is labeled properly)

[fr, time] = lowpass_filt_firingRate(mea, 1000);  % show the recruitment wave downsampled to ~100 Hz
inds_50_60 = logical( (time > 50) .* (time < 60) );  % ... only from 50-60 s
fr = fr(inds_50_60, :);  % ... resize to only this period
time = time(inds_50_60);  % ... same

[~, mea] = filter_mea(mea, [], {'lfp'; 'mua'});  % filter the raw data
event_inds = mua_events(mea);  % compute event times

inds_lfp = interp1(downsample(Time, mea.skipfactor), ...  % Select subset of lfp corresponding to 50-60 s
	1:length(mea.lfp), time, 'nearest');  
	
ds = abs(diff(time(1:2))) * mea.SamplingRate;  % compute the binwidth
bincounts = histcounts(event_inds, (1:ds:prod(size(mea.mua))) + ds/2);  % ... and the counts per bin
fr_high = reshape(bincounts, [], size(mea.mua, 2));  % Convert to array

close all;
framerate = 100;
electrode_vid({fr; fr_high(inds_50_60, :); mea.lfp(inds_lfp, :)}, ...
	P(:, 1), P(:, 2), time, time(end), framerate, ...
	mea.Name, true)
