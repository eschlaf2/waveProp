% Creates a video of the electrode array for the matfile variable mea,
% which should be in the workspace already.
% Example:
%     mea = matfile('c5/c5_Seizure1_Neuroport.mat');
%     electrode_vid_full;

[fr, time] = lowpass_filt_firingRate(mea);
try
	lfp = mea.lfp;
catch 
	lfp = filter_mea(mea, [], {'lfp'});
	lfp = lfp.lfp;
end
Time = mea.Time;
Time = Time();
te = Time(end) - mea.Padding(1, 2);
inds = arrayfun(@(t) find(Time == t, 1), time);
if isprop(mea, 'skipfactor')
	skipfactor = mea.skipfactor;
else
	skipfactor = 1;
end
	
fr_high = zeros(size(mea.mua));
event_inds = mea.event_inds;
t0 = find(Time == 0, 1);
event_inds(event_inds < t0) = [];
fr_high(event_inds) = 1;
ds = inds(2) - inds(1);
fr_high = smoothdata(fr_high, 1, 'movmean', ds) * ds;
close all;
electrode_vid({fr; fr_high(inds, :); lfp(round(inds/skipfactor) + 1, :)}, ...
	mea.X, mea.Y, time, time(end), [], ...
	mea.Name, true)