% Creates a video of the electrode array for the matfile variable mea,
% which should be in the workspace already. Limits video to seconds a to b
% Example:
%     mea = matfile('c5/c5_Seizure1_Neuroport.mat');
%     electrode_vid_full;

patpath = genpath(pat);
addpath(patpath);

fname = sprintf('%s_Seizure%d_Neuroport', pat, seizure);
if ~exist([fname '.mat'], 'file')
	if exist([fname '_10_10.mat'], 'file')
		fname = [fname '_10_10.mat'];
		fprintf('Using %s\n', fname);
	elseif exist(sprintf('%s_Seizure%d.mat', pat, seizure), 'file')
		disp('Creating epoch file ...')
		create_epoch(pat, seizure, 'padding', [10 10]);
		fprintf('Using %s\n', fname);
		fname = [fname '_10_10.mat'];
	end
end

mea = load(fname);

% mea = load('c5/c5_Seizure1_Neuroport_10_10.mat');  % load data
mea = exclude_channels(mea);  % remove non-participating channels
P = mea.Position;  % get electrode coordinates
P(mea.BadChannels, :) = [];  % ... and remove bad channels
Time = mea.Time;  % Get sampling times
Time = Time();  % ... (and convert from function to array if necessary)
te = Time(end) - mea.Padding(1, 2);  % Find when seizure ends (so that video is labeled properly)

[fr, time] = lowpass_filt_firingRate(mea, 1000);  % show the recruitment wave downsampled to ~100 Hz
inds_a_b = logical( (time > a) .* (time < b) );  % ... only from a-b s
fr = fr(inds_a_b, :);  % ... resize to only this period
time = time(inds_a_b);  % ... same

[~, mea] = filter_mea(mea, [], {'lfp'; 'mua'});  % filter the raw data
event_inds = mua_events(mea);  % compute event times

[nT, nCh] = size(mea.mua);
ds = round(abs(diff(time(1:2))) * mea.SamplingRate);  % compute the binwidth
fr_high = false([nT, nCh]);
fr_high(event_inds) = true;
fr_high = fr_high(1:(floor(nT / ds) * ds), :);  % trim
fr_high = squeeze(sum(reshape(fr_high, ds, [], nCh)));

inds_lfp = interp1(downsample(Time, mea.skipfactor), ...  % Select subset of lfp corresponding to 50-60 s
	1:length(mea.lfp), time, 'nearest');  
inds_fr = interp1(Time((1:length(fr_high)) * ds), ...
	1:length(fr_high), time, 'nearest');

close all;
framerate = 100;
electrode_vid({fr; fr_high(inds_fr, :); mea.lfp(inds_lfp, :)}, ...
	P(:, 1), P(:, 2), time, time(end), framerate, ...
	sprintf('%s_%d_%d', mea.Name, a, b), true)
