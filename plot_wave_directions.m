function output = plot_wave_directions(mea, waveFit, sig)
% Make a summary plot of wave directions 

% Figure subplot dimensions
rows = 2;
cols = 10; 

%% Parse input and set defaults
if ~exist('sig', 'var')
	sig = 0.05;
end

plotTitle = strrep(waveFit.Name, '_', ' ');
Z = waveFit.Z;  
V = waveFit.V;
fr = mea.firingRate;
meanFr = mean(fr, 2, 'omitnan');
Time = mea.Time;
Time = Time();
numWaves = numel(Z);
computeTimes = waveFit.computeTimes / 1e3;

%% Plot the firing rate and wave velocity at each discharge
output(1) = figure(); clf; fullwidth(true)
p(1) = subplot(rows, cols, 1:7);  % Left plot
p(2) = subplot(rows, cols, 9:10); % right plot (compass)
title(p(1), plotTitle)
xlabel(p(1), 'Time (s)')

% Left plot, right axis
% Put wave direction labels on the right y-axis	
yyaxis(p(1), 'right');  
ylim(p(1), [-pi, pi]);  % set limits
hold(p(1), 'on');
yticks(p(1), [-pi, 0, pi]);  % Put ticks at pi rad
ytks = p(1).YTick;  % store ticks (transform them to left axis)

% Relabel the right axis using arrows
% ylab = get(p(1), 'yticklabels');
labs = {'\leftarrow'; '\rightarrow'; '\leftarrow'};
p(1).YTickLabel = labs;
% ylabel(p(1), 'Wave direction (\pi rad)')

% Transform Z and ytks so that you can use the left axis
% ... (this way you get a colorbar)
x = [1 -pi; 1 pi];
y = [min(meanFr); max(meanFr)];
b = x\y;

% Left plot, left axis
% Create grid lines corresponding to the right axis
yyaxis(p(1), 'left')
i = 1;
tempc = lines(2);  % grid line colors correspond to angle
for y = ytks * pi
	i = mod(i + 1, 2);
	plot(p(1), [Time(1); Time(end)],  b' * [1; y] * [1; 1], ... 
		':', 'color', tempc(1 + i, :))
	hold on
end
xlim(p(1), [Time(1), Time(end)])
ylim(p(1), x * b);

% Plot the mean firing rate
plot(p(1), Time, meanFr, 'color', .5*[1 1 1]); 
ylabel(p(1), 'Mean firing rate (spikes/s)')
p(1).XGrid = 'on';
hold(p(1), 'on');

% Right plot
% Plot the velocities in a compass
hold(p(2), 'off')
compass(p(2), V(1, waveFit.p < sig), V(2, waveFit.p < sig)); 
hold(p(2), 'on');
colormap(p(2), 'cool');


cmap = cool(numWaves);  % compass color corresponds to time
cmapDir = hsv;  % scatter color corresponds to direction
colormap(p(1), cmapDir);
p(1).CLim = [-pi pi];
	
for i = 1:numWaves  % Overlay colored data points each discharge
	if waveFit.p(i) > sig
		continue
	end
	yyaxis(p(1), 'left'); % right axis
	scatter(p(1), computeTimes(i), b(1) + b(2) * Z(i), 40, waveFit.Z(i), 'filled');
	h = compass(p(2), waveFit.V(1, i), waveFit.V(2, i));
	h.Color = cmap(i, :);
	h.LineWidth = 2;

% 		drawnow()
% 		pause(1e-2)
end
hold(p(1), 'off')
hold(p(2), 'off')

% Add colorbars
p(2).CLim = [computeTimes(1) computeTimes(end)];
colorbar(p(2), 'southoutside', 'Ticks', [computeTimes(1) computeTimes(end)], ...
	'TickLabels', {'start'; 'end'})

%% Plot direction and speed as a function of time
p(3) = subplot(rows, cols, 11:17);
speed = sqrt(sum(V.^2));
stem(p(3), computeTimes, log(speed), 'linewidth', 2);
hold on
scatter(p(3), computeTimes, log(speed), 40, waveFit.Z, 'filled');
ylabel('Wave speed (log)')
colormap(p(3), cmapDir);
xlabel(p(3), 'Time (s)')
grid on

%% Create histograms of directions in first and last n seconds
n = 20;
firstInds = logical((computeTimes >= 0) .* (computeTimes <= n));
te = Time(end) - mea.Padding(2);
lastInds = logical((computeTimes >= te - n) .* (computeTimes <= te));
p(4) = subplot(rows, cols, 19:20);
hrose = rose(waveFit.Z(firstInds)); 
hrose.Color = tempc(1, :); 
hrose.LineWidth = 2; 
title(sprintf('Discharges in the first and last %d seconds\n', n))

hold(p(4), 'on');
hroseR = rose(waveFit.Z(lastInds)); 
hroseR.Color = tempc(end, :); 
hroseR.LineWidth = 2; 
axis(p(4), 'tight')

% Fix positioning
drawnow();
p(3).Position([1 3 4]) = p(1).Position([1 3 4]);
p(3).XLim = p(1).XLim;
c = colorbar(p(3), 'northoutside');
c.Label.String = 'Wave direction (rad)';

output = {output; p};