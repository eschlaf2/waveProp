function output = plot_wave_directions(mea, method, sig)
% Make a summary plot of wave directions computed according to method (NYC
% or BOS).

%% Parse input and set defaults
if ~exist('sig', 'var')
	sig = 0.05;
end
if ~exist('method', 'var') || isempty(method)  % Default to NYC method
	method = 'nyc';
end

switch lower(method)
	case 'bos'
		wave_fit = mea.wave_fit_bos;
		wave_fit.V = wave_fit.V';
		waveTimes = wave_fit.wave_times;
		figbase = 0;
	case 'nyc'
		wave_fit = mea.wave_fit_nyc;
		waveTimes = mea.waveTimes;
		figbase = 10;
end

plotTitle = [strrep(mea.Name, '_', ' ') ' ' upper(method) ' method'];
Z = wave_fit.Z;  
V = wave_fit.V;
fr = mea.firingRate;
meanFr = mean(fr, 2, 'omitnan');
Time = mea.Time;
TimeMs = Time() * 1e3;
numWaves = numel(Z);

%% Plot the firing rate and wave velocity at each discharge
output.main = figure(figbase + 2); clf; fullwidth()
p1 = subplot(1, 10, 1:7);  % Left plot
p2 = subplot(1, 10, 9:10); % right plot (compass)
% 	p2.NextPlot = 'replacechildren';
title(p1, plotTitle)
xlabel(p1, 'Time (ms)')

% Left plot, right axis
% Put wave direction labels on the right y-axis	
yyaxis(p1, 'right');  
ylim(p1, [min(Z) / pi, max(Z) / pi]);  % set limits
hold(p1, 'on');
yticks(p1, (floor(min(Z) / pi) : ceil(max(Z) / pi)));  % Put ticks at pi rad
ytks = p1.YTick;  % store ticks (transform them to left axis)

% Relabel the right axis using arrows
ylab = get(p1, 'yticklabels');
labs = {'\rightarrow', '\leftarrow'};
p1.YTickLabel = cellfun(@(x) labs{mod(round(10 * str2double(x)) / 10, 2) + 1}, ylab, 'uni', 0);
ylabel(p1, 'Wave direction (\pi rad)')

% Transform Zu and ytks so that you can use the left axis
% ... (this way you get a colorbar)
x = [1 min(Z); 1 max(Z)];
y = [min(meanFr); max(meanFr)];
b = x\y;

% Left plot, left axis
% Create grid lines corresponding to the right axis
yyaxis(p1, 'left')
i = 1;
tempc = lines(2);  % grid line colors correspond to angle
for y = ytks * pi
	i = mod(i + 1, 2);
	plot(p1, [TimeMs(1); TimeMs(end)],  b' * [1; y] * [1; 1], ... 
		':', 'color', tempc(1 + i, :))
	hold on
end
xlim(p1, [TimeMs(1), TimeMs(end)])
ylim(p1, x * b);

% Plot the mean firing rate
plot(p1, TimeMs, meanFr, 'color', .5*[1 1 1]); 
ylabel(p1, 'Mean firing rate (spikes/s)')
p1.XGrid = 'on';
hold(p1, 'on');

% Right plot
% Plot the velocities in a compass
hold(p2, 'off')
compass(p2, V(1, wave_fit.p < sig), V(2, wave_fit.p < sig)); 
hold(p2, 'on');
colormap(p2, 'cool');


cmap = cool(numWaves);  % compass color corresponds to time
cmapDir = hsv;  % scatter color corresponds to direction

for i = 1:numWaves  % Overlay colored data points each discharge
	if wave_fit.p(i) > sig
		continue
	end
	yyaxis(p1, 'left'); % right axis
	scatter(p1, waveTimes(i), b(1) + b(2) * Z(i), 30, wave_fit.Z(i), 'filled');
	colormap(p1, cmapDir);
	h = compass(p2, wave_fit.V(1, i), wave_fit.V(2, i));
	h.Color = cmap(i, :);
	h.LineWidth = 2;

% 		drawnow()
% 		pause(1e-2)
end
hold(p1, 'off')
hold(p2, 'off')

% Add colorbars
colorbar(p1)
p2.CLim = [waveTimes(1) waveTimes(end)];
colorbar(p2, 'southoutside', 'Ticks', [waveTimes(1) waveTimes(end)], ...
	'TickLabels', {'start'; 'end'})

% Create histograms of first and last n discharges
% if feature('ShowFigureWindows')
	n = 20;
	first_wave = find(mea.waveTimes > 0, 1);
	output.first = figure(figbase + 5);  % Plot a histogram of the first n discharges
	hrose = rose(wave_fit.Z(first_wave : first_wave + n)); 
	hrose.Color = tempc(1, :); 
	hrose.LineWidth = 2; 
	title('Direction during first 20 discharges')

	output.last = figure(figbase + 6);  % ... and the last twenty discharges
	hroseR = rose(wave_fit.Z(end - n : end)); 
	hroseR.Color = tempc(end, :); 
	hroseR.LineWidth = 2; 
	title('Direction during last 20 discharges')
% end