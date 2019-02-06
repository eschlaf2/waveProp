function [mea] = wave_prop(mea, PLOT)

%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency

if ~exist('PLOT', 'var')
	PLOT = true;
end

sig = .25;
halfWin = 50;  % half window around discharge event (ms)
plotTitle = strrep(mea.Name, '_', ' ');

%% Find discharge times

if ~isfield(mea, 'waveTimes')
	[~, waveTimes] = findpeaks(mean(mea.firingRate, 2), ...  % find peaks in mean firing rate
		mea.SamplingRate / 1e3, ...  % ... in ms 
		'minpeakprom', 10, ...  % ... use discrete peaks
		'minpeakdistance', 100);  % ... peaks should be at least 100 ms apart

	waveTimes = waveTimes - mea.Padding(1) * 1e3;  % Account for padding (in ms)
	mea.TimeMs = mea.Time * 1000;  % Convert times to ms

	% Sizing variables
	numCh = numel(mea.X);
	numWaves = numel(waveTimes);
	mea.waveTimes = waveTimes;
end

% Create an array of spike times
T = mea.TimeMs' .* mea.events;
T(T == 0) = nan;

	%% Estimate wave direction at each discharge time

if ~isfield(mea, 'fit')
	% Initialize arrays
	beta = nan(3, numWaves);  % fit parameters
	V = nan(2, numWaves);  % wave velocity (psuedo-inverse of beta)
	p = nan(1, numWaves);  % certainty

	for i = 1:numWaves  % estimate wave velocity for each discharge
		t = mea.waveTimes(i);
		inds = find((mea.TimeMs >= t - halfWin) .* (mea.TimeMs <= t + halfWin));
		events = T(inds, :)';
		events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
		for ch = 1:numCh
			temp = events{ch};
			temp(isnan(temp)) = [];
			events{ch} = temp;
		end
		[beta(:, i), V(:, i), p(i)] = ...
			SpatialLinearRegression(events, [mea.X, mea.Y], ...
			'switch_plot', 0, 'Lossfun','L2');
	end
	Z = angle(complex(V(1, :), V(2, :)));
	Zu = unwrap(Z);

	mea.fit.beta = beta;
	mea.fit.V = V;
	mea.fit.p = p;
	mea.fit.Z = Z;
	mea.fit.Zu = Zu;
end

%% Plot results
% Plot the mean firing rate along with the unwrapped wave angle at each
% discharge event.

if PLOT

	figure(2); clf; fullwidth()
	p1 = subplot(1, 10, 1:7);  % Left plot
	
	
	yyaxis(p1, 'right');
	ylim([min(mea.fit.Zu), max(mea.fit.Zu)]); hold on;
	yticks(floor(min(mea.fit.Zu) / pi) * pi : pi : ceil(max(mea.fit.Zu) / pi) * pi);
	
	for y = p1.YTick
		plot(p1, [mea.TimeMs(1); mea.TimeMs(end)],  y * [1; 1], ':', 'color', p1.YLabel.Color)
	end
	xlim([mea.TimeMs(1), mea.TimeMs(end)])
	
	plot(p1, mea.waveTimes, mea.fit.Z, 'k.')
	
	ylab = get(p1, 'yticklabels');
	labs = {'\rightarrow', '\leftarrow'};
	p1.YTickLabel = cellfun(@(x) labs{mod(round(10 * str2double(x) / pi) / 10, 2) + 1}, ylab, 'uni', 0);
	
	ylabel('Wave direction (\pi rad)')
	title(plotTitle)
	xlabel('Time (ms)')
	
	yyaxis left  % Plot the mean firing rate
	plot(p1, mea.TimeMs, mean(mea.firingRate, 2), 'color', .5*[1 1 1]); 
	ylabel('Mean firing rate (spikes/s)')

	p2 = subplot(1, 10, 9:10); 
	compass(mea.fit.V(1, mea.fit.p < sig), mea.fit.V(2, mea.fit.p < sig)); hold on;

	cmap = jet(numWaves); 
	for i = 2:numWaves
		if mea.fit.p(i) > sig
			continue
		end
% 		get(p1);  % left plot
		yyaxis(p1, 'right'); % right axis
		plot(p1, mea.waveTimes(i), mea.fit.Zu(i), 'ko', 'MarkerFaceColor', cmap(i, :));

% 		get(p2); 
		h = compass(p2, mea.fit.V(1, i), mea.fit.V(2, i)); hold on
		h.Color = cmap(i-1, :);
		h.LineWidth = 2;

		drawnow()
		pause(1e-2)
	end
	hold off
	
end



