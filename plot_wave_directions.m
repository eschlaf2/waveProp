function output = plot_wave_directions(mea, waveFit, sig, ax)
% Make a summary plot of wave directions 

%% Parse input and set defaults
if ~exist('sig', 'var'), sig = 0.05; end
if ~exist('ax', 'var'), ax = axes(); end

plotTitle = strrep(waveFit.Name, '_', ' ');
Z = waveFit.Z;  
try fr = mea.firingRate; catch; [fr, mea] = mua_firing_rate(mea); end
assignin('base', 'mea', mea);
meanFr = mean(fr, 2, 'omitnan');
Time = mea.Time;
Time = Time();
numWaves = numel(Z);
computeTimes = waveFit.computeTimes / 1e3;

%% Plot the firing rate and wave velocity at each discharge
    
title(ax, plotTitle)
xlabel(ax, 'Time (s)')

% Left plot, right axis
% Put wave direction labels on the right y-axis	
yyaxis(ax, 'right');  
ylim(ax, [-pi, pi]);  % set limits
hold(ax, 'on');
yticks(ax, [-pi, 0, pi]);  % Put ticks at pi rad
ytks = ax.YTick;  % store ticks (transform them to left axis)

% Relabel the right axis using arrows
% ylab = get(p, 'yticklabels');
labs = {'\leftarrow'; '\rightarrow'; '\leftarrow'};
ax.YTickLabel = labs;
% ylabel(p, 'Wave direction (\pi rad)')

% Transform Z and ytks so that you can use the left axis
% ... (this way you get a colorbar)
x = [1 -pi; 1 pi];
y = [min(meanFr); max(meanFr)];
b = x\y;

% Left plot, left axis
% Create grid lines corresponding to the right axis
yyaxis(ax, 'left')
i = 1;
tempc = lines(2);  % grid line colors correspond to angle
for y = ytks * pi
	i = mod(i + 1, 2);
	plot(ax, [Time(1); Time(end)],  b' * [1; y] * [1; 1], ... 
		':', 'color', tempc(1 + i, :))
	hold on
end
xlim(ax, [Time(1), Time(end)])
ylim(ax, x * b);

% Plot the mean firing rate
plot(ax, Time, meanFr, 'color', .5*[1 1 1]); 
ylabel(ax, 'Mean firing rate (spikes/s)')
ax.XGrid = 'on';
hold(ax, 'on');

cmapDir = hsv;  % scatter color corresponds to direction
colormap(ax, cmapDir);
ax.CLim = [-pi pi];
% scatter(ax, computeTimes, b(1) + b(2) * Z, 40, 'filled', 'MarkerFaceAlpha', 1);
	
for i = 1:numWaves  % Overlay colored data points each discharge
	if waveFit.p(i) >= sig
		continue
	end
	yyaxis(ax, 'left'); % right axis
	scatter(ax, computeTimes(i), b(1) + b(2) * Z(i), 40, waveFit.Z(i), 'filled');

% 		drawnow()
% 		pause(1e-2)
end
hold(ax, 'off')


%% return
if nargout > 0
	output.ax = ax;
	output.meanFr = meanFr;
	output.b = b;
end