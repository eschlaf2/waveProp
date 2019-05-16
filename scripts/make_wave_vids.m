% mea = ...
% [fit, mea] = wave_prop(mea, ...)

% fit = fit_dev;
time = downsample(mea.Time(), mea.skipfactor) * 1e3;
position = mea.Position;
try
	position(mea.BadChannels, :) = [];
catch ME
end

numwaves = length(fit.computeTimes);
halfwin = 50;  % ms
inds = @(t) logical((time >= (t - halfwin)) .* (time <= (t + halfwin)));

for ii = 1:numwaves
	v = VideoWriter(sprintf('%s_wave%03d', mea.Name, ii));
	v.FrameRate = 100;
	open(v);
	h = figure(1); clf; ax = axes();
	title(ax, sprintf('Wave %d', ii));
	set(h, 'nextplot', 'replacechildren')
	plot_window(mea.lfp(inds(fit.computeTimes(ii)), :), position, ax, v);
	close(v)
	fprintf('%d/%d done\n', ii, numwaves)
end