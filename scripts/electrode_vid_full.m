% Creates a video of the electrode array for the matfile variable mea,
% which should be in the workspace already.
% Example:
%     mea = matfile('c5/c5_Seizure1_Neuroport.mat');
%     electrode_vid_full;

mea = exclude_channels(mea);  % Exclude channels based on activity

dsRate = 100;  % Set the desired downsampling rate for the video (frameRate is half of this... below)
[fr, time] = lowpass_filt_firingRate(mea);  % create the data for the left plot
[~, mea] = filter_mea(mea, {'lfp'; 'mua'});  % filter the raw data

Time = mea.Time;  % Get relevant time points
Time = Time();
te = Time(end) - mea.Padding(1, 2);
t0 = find(Time == 0, 1);
inds = arrayfun(@(t) find(Time == t, 1), time);

P = mea.Position;  % Get electrode coordinates
P(mea.BadChannels, :) = [];  % remove excluded channels
	
% Get MUA activity
fr_high = zeros(size(mua));
event_inds = mua_events(mea);

fr_high(event_inds) = 1;
ds = inds(2) - inds(1);
fr_high = smoothdata(fr_high, 1, 'movmean', ds) * ds;
close all;
electrode_vid({fr; fr_high(inds, :); lfp(round(inds/skipfactor) + 1, :)}, ...
	P(:, 1), P(:, 2), time, time(end), dsRate / 2, ...
	mea.Name, true)