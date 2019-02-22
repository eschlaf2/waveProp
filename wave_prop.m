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

try 
	waveTimes = mea.waveTimes;
catch ME
	[~, waveTimes] = findpeaks(mean(mea.firingRate, 2), ...  % find peaks in mean firing rate
		mea.SamplingRate / 1e3, ...  % ... in ms 
		'minpeakprom', 10, ...  % ... use discrete peaks
		'minpeakdistance', 100);  % ... peaks should be at least 100 ms apart

	padding = mea.Padding;
	waveTimes = waveTimes - padding(1) * 1e3;  % Account for padding (in ms)
	TimeMs = mea.Time * 1000;  % Convert times to ms

	% Sizing variables
	numCh = numel(mea.X);
	numWaves = numel(waveTimes);
	
	% Save wave times
	mea.waveTimes = waveTimes;
end

% Create an array of spike times
T = nan(size(mea.mua), 'single');
T(mea.event_inds) = 1;
T = TimeMs' .* T;
% T(T == 0) = nan;

%% Estimate wave direction at each discharge time

try 
	wave_fit = mea.fit;
catch
	% Initialize arrays
	beta = nan(3, numWaves);  % fit parameters
	V = nan(2, numWaves);  % wave velocity (psuedo-inverse of beta)
	p = nan(1, numWaves);  % certainty
	position = [mea.X, mea.Y];

	for i = 1:numWaves  % estimate wave velocity for each discharge
		t = waveTimes(i);
		inds = find((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
		events = T(inds, :)';
		events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
		for ch = 1:numCh
			temp = events{ch};
			temp(isnan(temp)) = [];
			events{ch} = temp;
		end
		[beta(:, i), V(:, i), p(i)] = ...
			SpatialLinearRegression(events, position, ...
			'switch_plot', 0, 'Lossfun','L2');
	end
	Z = angle(complex(V(1, :), V(2, :)));
	Zu = unwrap(Z);

	wave_fit.beta = beta;
	wave_fit.V = V;
	wave_fit.p = p;
	wave_fit.Z = Z;
	wave_fit.Zu = Zu;
	mea.fit = wave_fit;
end

%% Plot results
% Plot the mean firing rate along with the unwrapped wave angle at each
% discharge event.
% set(2,'DefaultAxesColorOrder', get(groot, 'factoryAxesColorOrder'));

if PLOT
	
	% Plot the firing rate and wave velocity at each discharge
	figure(2); clf; fullwidth()
	p1 = subplot(1, 10, 1:7);  % Left plot
	p2 = subplot(1, 10, 9:10); % right plot (compass)
	title(plotTitle)
	xlabel('Time (ms)')
	
	% Left plot, right axis
	% Put wave direction labels on the right y-axis	
	yyaxis(p1, 'right');  
	ylim([min(wave_fit.Zu) / pi, max(wave_fit.Zu) / pi]);  % set limits
	hold on;
	yticks(p1, (floor(min(wave_fit.Zu) / pi) : ceil(max(wave_fit.Zu) / pi)));  % Put ticks at pi rad
	ytks = p1.YTick;  % store ticks (transform them to left axis)
	
	% Relabel the right axis using arrows
	ylab = get(p1, 'yticklabels');
	labs = {'\leftarrow', '\rightarrow'};
	p1.YTickLabel = cellfun(@(x) labs{mod(round(10 * str2double(x)) / 10, 2) + 1}, ylab, 'uni', 0);
	ylabel('Wave direction (\pi rad)')

    % Transform Zu and ytks so that you can use the left axis
	% ... (this way you get a colorbar)
	fr = mean(mea.firingRate, 2);
	x = [1 min(wave_fit.Zu); 1 max(wave_fit.Zu)];
	y = [min(fr); max(fr)];
	b = x\y;
	
	% Left plot, left axis
	% Create grid lines corresponding to the right axis
	yyaxis(p1, 'left')
	i = 1;
	tempc = hsv(2);  % grid line colors correspond to angle
	for y = ytks * pi
		i = mod(i + 1, 2);
		plot(p1, [TimeMs(1); TimeMs(end)],  b' * [1; y] * [1; 1], ... 
			':', 'color', tempc(1 + i, :))
		hold on
	end
	xlim([TimeMs(1), TimeMs(end)])
	ylim(x * b);
	
	% Plot the mean firing rate
	plot(p1, TimeMs, fr, 'color', .5*[1 1 1]); 
	ylabel('Mean firing rate (spikes/s)')
	
	% Right plot
	% Plot the velocities in a compass
	compass(p2, wave_fit.V(1, wave_fit.p < sig), wave_fit.V(2, wave_fit.p < sig)); hold on;
	colormap('cool');

	
	cmap = cool(numWaves);  % compass color corresponds to time
	cmapDir = hsv;  % scatter color corresponds to direction
	
	for i = 2:numWaves  % Overlay colored data points each discharge
		if wave_fit.p(i) > sig
			continue
		end
		yyaxis(p1, 'left'); % right axis
		scatter(p1, waveTimes(i), b(1) + b(2) * wave_fit.Zu(i), 30, wave_fit.Z(i), 'filled');
		colormap(p1, cmapDir);
		h = compass(p2, wave_fit.V(1, i), wave_fit.V(2, i)); hold on
		h.Color = cmap(i-1, :);
		h.LineWidth = 2;

% 		drawnow()
% 		pause(1e-2)
	end
	hold off
	
	% Add colorbars
	colorbar(p1)
	p2.CLim = [waveTimes(1) waveTimes(end)];
	colorbar(p2, 'southoutside', 'Ticks', [waveTimes(1) waveTimes(end)], ...
		'TickLabels', {'start'; 'end'})
	
	% Create histograms of first and last n discharges
	if feature('ShowFigureWindows')
		n = 20;
		figure(5);  % Plot a histogram of the first twenty discharges
		hrose = rose(wave_fit.Z(5 : 5 + n)); 
		hrose.Color = temp(1, :); 
		hrose.LineWidth = 2; 
		title('Direction during first 20 discharges')
		
		figure(6);  % ... and the last twenty discharges
		hroseR = rose(wave_fit.Z(end - n : end)); 
		hroseR.Color = temp(end, :); 
		hroseR.LineWidth = 2; 
		title('Direction during last 20 discharges')
	end
end



