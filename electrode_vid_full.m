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
	skipfactor = lfp.skipfactor;
	lfp = lfp.lfp;
end
try
	mua = mea.mua;
catch 
	mua = filter_mea(mea, [], {'mua'});
	mua = mua.mua;
end

Time = mea.Time;
Time = Time();
te = Time(end) - mea.Padding(1, 2);
inds = arrayfun(@(t) find(Time == t, 1), time);
	
fr_high = zeros(size(mua));
if any(strcmpi(fieldnames(mea), 'event_inds'))
	event_inds = mea.event_inds;
else
	event_inds = mua_events(mea);
end
t0 = find(Time == 0, 1);
event_inds(event_inds < t0) = [];
fr_high(event_inds) = 1;
ds = inds(2) - inds(1);
fr_high = smoothdata(fr_high, 1, 'movmean', ds) * ds;
close all;
electrode_vid({fr; fr_high(inds, :); lfp(round(inds/skipfactor) + 1, :)}, ...
	mea.X, mea.Y, time, time(end), [], ...
	mea.Name, true)