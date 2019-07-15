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
try fr = mea.firingRate; catch; [fr, mea] = mua_firing_rate(mea); end
assignin('base', 'mea', mea);
meanFr = mean(fr, 2, 'omitnan');
Time = mea.Time;
Time = Time();
numWaves = numel(Z);
computeTimes = waveFit.computeTimes / 1e3;

%% Plot the firing rate and wave velocity at each discharge

if 1
    
output(1) = figure(); clf; fullwidth(true)
p(1) = subplot(rows, cols, 1:7);  % Left plot
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

% 		drawnow()
% 		pause(1e-2)
end
hold(p(1), 'off')



end

%% Compass plot

if 0
    
p(2) = subplot(rows, cols, 9:10); % right plot (compass)


% Right plot
% Plot the velocities in a compass
hold(p(2), 'off')
compass(p(2), V(1, waveFit.p < sig), V(2, waveFit.p < sig)); 
hold(p(2), 'on');
colormap(p(2), 'cool');


for i = 1:numWaves  % Overlay colored data points each discharge
	if waveFit.p(i) > sig
		continue
	end
	h = compass(p(2), waveFit.V(1, i), waveFit.V(2, i));
	h.Color = cmap(i, :);
	h.LineWidth = 2;

% 		drawnow()
% 		pause(1e-2)
end

hold(p(2), 'off')

% Add colorbars
p(2).CLim = [computeTimes(1) computeTimes(end)];
colorbar(p(2), 'southoutside', 'Ticks', [computeTimes(1) computeTimes(end)], ...
	'TickLabels', {'start'; 'end'})

end
%% Plot direction and speed as a function of time
if 0
    
p(3) = subplot(rows, cols, 11:17);
speed = ((log(sqrt(sum(V.^2)))));
zscore = @(x) (x - mean(x, 'omitnan')) / std(x, 'omitnan');
stem(p(3), computeTimes, zscore(speed), 'linewidth', 1);
hold on
scatter(p(3), computeTimes, zscore(speed), 40, waveFit.Z, 'filled');
ylabel('Wave speed (log)')
colormap(p(3), cmapDir);
p(3).CLim = [-pi pi];
ylim(p(3), [-4, 4]);
xlabel(p(3), 'Time (s)')
grid on

end

%% Create histograms of directions in first and last n seconds

if 0
    
n = 20;
firstInds = logical((computeTimes >= 0) .* (computeTimes <= n));
te = Time(end) - mea.Padding(2);
lastInds = logical((computeTimes >= te - n) .* (computeTimes <= te));
edges = linspace(-pi, pi, 13);
p(4) = subplot(rows, cols, 19:20, 'replace');
p(4) = polarhistogram(waveFit.Z(logical(firstInds + lastInds)), edges, ...
	'FaceColor', 'none', 'LineStyle', 'none'); 
hold on;
r1 = polarhistogram(waveFit.Z(firstInds), edges); 
r1.FaceColor = cmap(1, :); 
r1.LineWidth = 2; 

r2 = polarhistogram(waveFit.Z(lastInds), edges); 
r2.FaceColor = cmap(end, :); 
r2.LineWidth = 2; 
% axis(p(4), 'tight')
title(sprintf('Discharges in the first and last %d seconds\n', n))


% Fix positioning
drawnow();
p(3).Position([1 3 4]) = p(1).Position([1 3 4]);
p(3).XLim = p(1).XLim;
c = colorbar(p(3), 'northoutside');
c.Label.String = 'Wave direction (rad)';

end


%% return
output = {output; p};