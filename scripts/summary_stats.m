% Plot the circular mean and standard deviation comparing chosen metrics
% for wave direction. Highlight individual patients and seizures.

% Load the table with comparison data
load('direction_stats.mat');  % This should have the stats and a variable with the metrics
pairs = [2 3];  % Ignore very short delay windows for now

% Extract variables from table
theta = stats.theta;
m1 = stats.m1;
cconf = stats.conf;
% cconf = real(cellfun(@(Z) circ_confmean(angle(Z(isfinite(Z)))), stats.dZ));
nanconf = isnan(cconf);
whichpair = stats.whichpair;
[S, C, sigS, sigC, cvar] = deal(nan(size(stats, 1)));
for ii = 1:size(stats, 1)
	S(ii) = mean(imag(stats.dZ{ii}), 'omitnan');
	C(ii) = mean(real(stats.dZ{ii}), 'omitnan');
	sigS(ii) = std(imag(stats.dZ{ii}), 'omitnan') / sqrt(stats.N(ii));
	sigC(ii) = std(real(stats.dZ{ii}), 'omitnan') / sqrt(stats.N(ii));
% 	r = sqrt(sigS(ii)^2 + sigC(ii)^2);
% 	r = max(sigS(ii), sigC(ii));
	Y = S(ii) + 2*sigS(ii) * [1 -1];
	X = C(ii) + 2*sigC(ii) * [1 -1];
% 	cvar(ii) = atan(2*r/stats.R(ii));
	[XX, YY] = meshgrid(linspace(X(1), X(2), 100), linspace(Y(1), Y(2), 100));	
	cvar(ii) = range(unwrap(sort(angle(complex(XX(:), YY(:)))))) / 2;
end
%%
% cconf(nanconf) = 2*stats.sigma(nanconf);
cconf(nanconf) = cvar(nanconf);
% cconf(nanconf) = pi;
% ... including variables for patient and seizure
[patient, seizure] = deal(cell(size(stats, 1), 1));
for ii = 1:size(stats, 1)
	info = strsplit(stats.filename{ii}, '_');
	patient{ii} = info{1};
	seizure{ii} = [' ' info{2}(8:end)];
end

% Move SIM to the end
mask = strcmpi(patient, 'sim');
reorder = @(X) [X(~mask); X(mask)];
patient = reorder(patient);
seizure = reorder(seizure);
m1 = reorder(m1);
theta = reorder(theta);
cconf = reorder(cconf);
nanconf = reorder(nanconf);
cvar = reorder(cvar);
whichpair = reorder(whichpair);

% Add a nan row between patients
[~, u] = unique(patient);
for uu = sort(u(u ~= 1), 'descend')'
	patient = [patient(1:uu-1); {''}; patient(uu:end)];
	seizure = [seizure(1:uu-1); {''}; seizure(uu:end)];
	theta = [theta(1:uu-1); nan; theta(uu:end)];
	cconf = [cconf(1:uu-1); nan; cconf(uu:end)];
	nanconf = [nanconf(1:uu-1); false; nanconf(uu:end)];
	cvar = [cvar(1:uu-1); nan; cvar(uu:end)];
	whichpair = [whichpair(1:uu-1); 0; whichpair(uu:end)];
	m1 = [m1(1:uu-1); 0; m1(uu:end)];
end

% Compute lower and upper bounds (circular variance)
lowCI = theta - cconf;
hiCI = theta + cconf;
toolow = lowCI < -pi;
toohi = hiCI > pi;

% Clear the figure
fig = figure(1); clf
ax = axes('parent', fig);
offset = @(ii) 2.2 * pi * ii;

% Plot the means first so it's easy to label
for ii = 1:length(pairs)
	mask = whichpair==ii | whichpair == 0;
	plot(theta(mask) + offset(ii), 'LineWidth', 2); hold on
end
set(ax, 'ColorOrderIndex', 1);

% Plot the gray blocks, confidence bounds, and means again
cdata = lines(7);
cicolor = .5 * [1 1 1];
for ii = 1:length(pairs)
	mask = whichpair==ii | whichpair == 0;
	N = sum(mask);
	
% 	if mod(ii, 2)
		patch([0 N+1 N+1 0], [-pi -pi pi pi] + offset(ii), .8*[1 1 1], ...
			'linestyle', 'none', 'facealpha', .5, 'linewidth', 2); hold on; 
% 	end
	x = [1;1] * (1:sum(mask));
	y = [lowCI(mask)'; hiCI(mask)'];
	toohi = y(2, :) > pi;
	toolow = y(1, :) < -pi;
	cicolor = cdata(ii, :);
	
	plot(ax, x, max(min(y, pi), -pi) + offset(ii), 'color', cicolor); 
	plot(ax, x(:, toohi), ...
		[0*y(1, toohi) - pi; y(2, toohi) - 2*pi] + offset(ii), ...
		'color', cicolor);
	plot(ax, x(:, toolow), ...
		[0*y(2, toolow) + pi; y(1, toolow) + 2*pi] + offset(ii), ...
		'color', cicolor); 
	scatter(ax, 1:N, theta(mask) + offset(ii), 50, cdata(ii, :), 'filled');
	x = 1:N;
	y = theta(mask);
	scatter(ax, x(nanconf(mask)), y(nanconf(mask)) + offset(ii), 90, .4*[1 1 1], 'd', 'filled');
% 	h.Marker
end
hold off

% % Put ticks at the center of each group (i.e. dZ=0)
% yticks((offset(1) - pi:offset(1):offset(ii)) + 1*pi)
% yticklabels([])  % ... but don't label since it will be in the legend
% 
% % Label each seizure
% xticks(1:N)
% xticklabels(seizure(mask))

% set(gca, 'ticklength', [0;0], 'ygrid', 'on')
axis tight
% box(ax, 'on')
set(ax,'TickLength',[0 0], ...
	'XTick', (1:N), ...
	'XTickLabel', seizure(mask), ...
	'YGrid','on', ...
	'YTick', (offset(1) - pi:offset(1):offset(ii)) + 1*pi, ...
	'YTickLabel', []);
% axis tight

%% make legend
metricpairs = strrep(metricpairs, 'maxdescent', 'M');
metricpairs = strrep(metricpairs, 'events', 'E');
metricpairs = strrep(metricpairs, 'delays_T01_fband1_50', 'D');
metricpairs = strrep(metricpairs, 'delays_T0p2_fband0_50', 'Ds');
l = cell(length(pairs), 1);
for ii = 1:length(pairs)
	l{ii} = sprintf('%s - %s \\color{white}...', metricpairs{pairs(ii), 2},metricpairs{pairs(ii), 1}); 
end
lgd = legend(l, ...
	'location', 'northoutside', ...
	'Orientation', 'horiz', ...
	'Interpreter', 'tex', ...
	'FontName', 'pt sans caption', ...
	'FontSize', 16);
lgd.Box = 'off';
ax.FontName = 'pt sans caption';
% ax.Box = 'off';
% ax.GridLineStyle = '-';
% ax.BoxStyle = 'full'
% ax.Line

