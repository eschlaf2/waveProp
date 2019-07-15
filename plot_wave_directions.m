function output = plot_wave_directions(mea, waveFit, sig)
% Make a summary plot of wave directions 

%% Parse input and set defaults
if ~exist('sig', 'var')
	sig = 0.05;
end

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
    
p = axes();
title(p, plotTitle)
xlabel(p, 'Time (s)')

% Left plot, right axis
% Put wave direction labels on the right y-axis	
yyaxis(p, 'right');  
ylim(p, [-pi, pi]);  % set limits
hold(p, 'on');
yticks(p, [-pi, 0, pi]);  % Put ticks at pi rad
ytks = p.YTick;  % store ticks (transform them to left axis)

% Relabel the right axis using arrows
% ylab = get(p, 'yticklabels');
labs = {'\leftarrow'; '\rightarrow'; '\leftarrow'};
p.YTickLabel = labs;
% ylabel(p, 'Wave direction (\pi rad)')

% Transform Z and ytks so that you can use the left axis
% ... (this way you get a colorbar)
x = [1 -pi; 1 pi];
y = [min(meanFr); max(meanFr)];
b = x\y;

% Left plot, left axis
% Create grid lines corresponding to the right axis
yyaxis(p, 'left')
i = 1;
tempc = lines(2);  % grid line colors correspond to angle
for y = ytks * pi
	i = mod(i + 1, 2);
	plot(p, [Time(1); Time(end)],  b' * [1; y] * [1; 1], ... 
		':', 'color', tempc(1 + i, :))
	hold on
end
xlim(p, [Time(1), Time(end)])
ylim(p, x * b);

% Plot the mean firing rate
plot(p, Time, meanFr, 'color', .5*[1 1 1]); 
ylabel(p, 'Mean firing rate (spikes/s)')
p.XGrid = 'on';
hold(p, 'on');

cmapDir = hsv;  % scatter color corresponds to direction
colormap(p, cmapDir);
p.CLim = [-pi pi];
	
for i = 1:numWaves  % Overlay colored data points each discharge
	if waveFit.p(i) >= sig
		continue
	end
	yyaxis(p, 'left'); % right axis
	scatter(p, computeTimes(i), b(1) + b(2) * Z(i), 40, waveFit.Z(i), 'filled');

% 		drawnow()
% 		pause(1e-2)
end
hold(p, 'off')


%% return
output = p;