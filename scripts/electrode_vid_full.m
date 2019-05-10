% Creates a video of the electrode array for the matfile variable mea,
% which should be in the workspace already.
% Example:
%     mea = matfile('c5/c5_Seizure1_Neuroport.mat');
%     electrode_vid_full;

% mea = load('/Users/emilyschlafly/Desktop/temp/c5_Seizure2_Neuroport_10_10.mat');
% mea_exc = exclude_channels(mea);  % Exclude channels based on activity

dsRate = 1000;  % Set the desired downsampled collection rate for the video (frameRate is half of this... below)
[~, mea] = filter_mea(mea, {'lfp'; 'mua'});  % filter the raw data

% [fr, time] = lowpass_filt_firingRate(mea, dsRate);  % create the data for the left plot
% Time = mea.Time;  % Get relevant time points
% Time = Time();
% te = Time(end) - mea.Padding(1, 2);
% t0 = find(Time == 0, 1);
% inds = arrayfun(@(t) find(Time == t, 1), time);

skipfactor = round(mea.SamplingRate / dsRate);
skiplfp = round(skipfactor / mea.skipfactor);
time = downsample(mea.Time(), skipfactor);
te = time(end) - mea.Padding(1, 2);
t0 = find(time >= 0, 1);

pos = mea.Position;  % Get electrode coordinates
pos(mea.BadChannels, :) = [];  % remove excluded channels

try
	[~, mea_exc] = filter_mea(mea_exc, {'lfp'; 'mua'});
	P_exc = mea.Position;
	P_exc(mea_exc.BadChannels, :) = [];
	data{5} = downsample(mea_exc.lfp, skiplfp); P{5} = P_exc;
catch ME
	if ~strcmp(ME.message, "Undefined function or variable 'mea_exc'.")
		rethrow(ME)
	end
end

%% Remove PCs
[coeff, score, latent, ~, explained, ~] = pca(mea.lfp);
k = find(cumsum(explained) > 90, 2);
pc = score(:, k(2):end) * coeff(:, k(2):end)';
data{3} = downsample(pc, skiplfp);
P{3} = pos;

%% Whiten
[~, Z] = zca_whitening(mea.lfp(time < 0, :));
data{1} = downsample(mea.lfp * Z, skiplfp);
P{1} = pos;

%% Remove mean
data{2} = downsample(mea.lfp - mean(mea.lfp, 2), skiplfp);
P{2} = pos;

%%
% Get MUA activity
% fr_high = zeros(size(mea.mua));
% event_inds = mua_events(mea);
% 
% fr_high(event_inds) = 1;
% ds = inds(2) - inds(1);
% fr_high = smoothdata(fr_high, 1, 'movmean', ds) * ds;
% data = {fr_high(inds, :), mea.lfp(round(inds / mea.skipfactor) + 1, :)};
data{4} = downsample(mea.lfp, skiplfp);
P{4} = pos;

data = data([2, 3, 4]); P = P([2,3,4]);
close all;
electrode_vid(data, P, time, te, 100, [mea.Name '_flat'], true);
surf_vid(data, P, time, te, 100, [mea.Name '_surf'], true);

% electrode_vid({fr; fr_high(inds, :); mea.lfp(round(inds/mea.skipfactor) + 1, :)}, ...
% 	P(:, 1), P(:, 2), time, time(end), dsRate / 2, ...
% 	mea.Name, true)