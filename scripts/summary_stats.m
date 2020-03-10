% Plot the circular mean and standard deviation comparing chosen metrics
% for wave direction. Highlight individual patients and seizures.

% Load the table with comparison data
% load('direction_stats.mat');  % This should have the stats and a variable with the metrics
% pairs = [2 3];  % Ignore very short delay windows for now
if ~exist('pairs', 'var'); pairs = 1:size(metricpairs, 1); end

% Extract variables from table
theta = stats.theta;
m1 = stats.m1;
cconf = stats.conf;
nanconf = isnan(cconf);
whichpair = stats.whichpair;

%%
cconf(nanconf) = pi;
% ... including variables for patient and seizure
[patient, seizure] = deal(cell(size(stats, 1), 1));
for ii = 1:size(stats, 1)
	patient{ii} = stats.name{ii}{1};
	seizure{ii} = stats.name{ii}{2};
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
% cvar = reorder(cvar);
whichpair = reorder(whichpair);

% Add a nan row between patients
[~, u] = unique(patient);
for uu = sort(u(u ~= 1), 'descend')'
	patient = [patient(1:uu-1); {''}; patient(uu:end)];
	seizure = [seizure(1:uu-1); {''}; seizure(uu:end)];
	theta = [theta(1:uu-1); nan; theta(uu:end)];
	cconf = [cconf(1:uu-1); nan; cconf(uu:end)];
	nanconf = [nanconf(1:uu-1); false; nanconf(uu:end)];
% 	cvar = [cvar(1:uu-1); nan; cvar(uu:end)];
	whichpair = [whichpair(1:uu-1); 0; whichpair(uu:end)];
	m1 = [m1(1:uu-1); 0; m1(uu:end)];
end

% Compute lower and upper bounds (circular variance)
lowCI = theta - cconf;
hiCI = theta + cconf;
toolow = lowCI < -pi;
toohi = hiCI > pi;

%% Make the figure
% Clear the figure
fig = figure(1); clf
nP = length(pairs);
ax = gobjects(nP, 1);

% Plot the means first so it's easy to label
cdata = lines(7);

% Rename metrics
metricpairs = rename_metrics(metricpairs);

for ii = 1:nP
	ax(ii) = subplot(nP, 1, ii);
	mask = whichpair==pairs(ii) | whichpair == 0;
	N = sum(mask);
	
% 	'XAxisLocation', 'origin', ...
	plot(theta(mask), 'LineWidth', 2, 'color', cdata(ii, :)); 
	
	set(ax(ii),'TickLength',[0 0], ...
		'XTick', (1:N), ...
		'XTickLabel', [], ...
		'YGrid','on', ...  
		'YTick', (-pi:pi/4:pi), ...
		'YTickLabel', [], ...
		'box', 'off', ...
		'nextplot', 'replacechildren', ...
		'xlim', [0 N+1], ...
		'ylim', [-pi pi]);
	hold(ax(ii), 'on')
end
set(ax(ii), 'XTickLabel', seizure(mask));

include_zero = lowCI <= 0 & hiCI >=0 & ~nanconf;
theta_small = abs(theta) < pi/4 & ~nanconf;

for ii = 1:nP
	mask = whichpair==pairs(ii) | whichpair == 0;
	N = sum(mask);
	
	x = [1;1] * (1:sum(mask));
	y = [lowCI(mask)'; hiCI(mask)'];
	toohi = y(2, :) > pi;
	toolow = y(1, :) < -pi;
	cicolor = cdata(ii, :);
	
	plot(ax(ii), x, y, 'color', cicolor); 
	plot(ax(ii), x(:, toohi), ...
		[0*y(1, toohi) - pi; y(2, toohi) - 2*pi], ...
		'color', cicolor);
	plot(ax(ii), x(:, toolow), ...
		[0*y(2, toolow) + pi; y(1, toolow) + 2*pi], ...
		'color', cicolor); 
	scatter(ax(ii), 1:N, theta(mask), 50, cdata(ii, :), 'filled');
	x = 1:N;
	y = theta(mask);
	scatter(ax(ii), ...
		x(nanconf(mask)), y(nanconf(mask)), 90, .4*[1 1 1], ...
		'd', 'filled');
	scatter(ax(ii), ...
		x(include_zero(mask)), y(include_zero(mask)), 90, [1 0 0]);  % indicate where difference includes zero with a circle
	plot(ax(ii), ...
		x(theta_small(mask)), y(theta_small(mask)), 'r*');  % indicate small differences with a *
	ylabel(ax(ii), sprintf('%s - %s \\color{white}...', ...
		metricpairs{pairs(ii), 2}, metricpairs{pairs(ii), 1} ...
		));
	hold(ax(ii), 'off');
	yline(ax(ii), 0);

end

str = annotation('textbox', ...
	'string', {'\circ := CI includes zero'; '* := \theta < \pi/4'}, ...
	'Position', [.01 .9 .05 .04], 'FitBoxToText', true);
