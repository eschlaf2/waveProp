data=load('testing/subsets/CUCX4_2_M.mat', '*stats');

f = fieldnames(data);
nF = numel(f);
dim = 1;  % dim=1: returns stats for each wave time; 
          % dim=2: returns states for each trial

% Convenience
mn =@(data) mean(data, dim, 'omitnan');
sd =@(data) std(data, [], dim, 'omitnan');
ci =@(data) 2 * sd(data) ./ sqrt(sum(isfinite(data), dim));

summary.dr = nan(nF, 2);
[summary.dspeed, summary.dtheta] = deal(nan(...  % Array of nans 
	nF, ...  % ... with first dimension for files,
	2, ...  % ... second dimension for mean and ci,
	size(data.(f{1}).dtheta, find([1 2] ~= dim)) ...  % third dimension for times or trials depending on dim
	));

for ii = 1:nF
	
	% Extract variables
	s = data.(f{ii});
	detections = s.detections;
	detection_rate = s.detection_rate;
	dtheta = s.dtheta;
	dspeed = s.dspeed;
	
	% Detection rates
	outliers = isoutlier(detection_rate);
	detection_rate(outliers) = nan;
	summary.dr(ii, 1) = mn(detection_rate);
	summary.dr(ii, 2) = ci(detection_rate);
	
	% Differences in speed
	outliers = isoutlier(dspeed, dim);
	dspeed(outliers) = nan;
	dspeed_mn = mn(dspeed);
	dspeed_mn(isoutlier(dspeed_mn)) = nan;
	summary.dspeed(ii, 1, :) = dspeed_mn;
	summary.dspeed(ii, 2, :) = ci(dspeed_mn);
	
	% Differences in direction
	outliers = isoutlier(abs(dtheta), dim);
	dtheta(outliers) = nan;
	summary.dtheta(ii, 1, :) = circ_mean(dtheta, [], dim, 'omitnan');
	confmean = circ_confmean(dtheta, [], [], [], dim, 'omitnan');
	if ~isreal(confmean), confmean = angle(confmean); end
	summary.dtheta(ii, 2, :) = confmean;
	
	
end

plot_detection_rates(summary, f);
plot_dthetaSpeed(summary.dtheta, f, '\Delta \theta');
plot_dthetaSpeed(summary.dspeed, f, '|\Delta V|');


%% Subroutines
function plot_dthetaSpeed(data, f, ttl)

% data = summary.dtheta;  % skip full data
% ttl = '\Delta\theta';

[nf, ~, nt] = size(data);

figure(); fullwidth(true);
nr = floor(sqrt(nf));
nc = ceil((nf) / nr);
x = 1:nt;

% Plot details
for ii = 1:nf
	subplot(nr, nc, ii)
	mn = squeeze(data(ii, 1, :));
	ci = squeeze(data(ii, 2, :));
	mn(isnan(ci)) = nan;
	
	plot(x .* [1; 1], mn' + [-1; 1] .* ci', 'r-', ...
		x, mn, '*');
	title(f{ii}(1:5));
	ylabel(ttl);
	
end

% Plot summary
mn_files = mean(data(:, 1, :), 3, 'omitnan');
std_files = std(data(:, 1, :), [], 3, 'omitnan');
n_obs_files = sum(isfinite(data(:, 1, :)), 3);
x = 1:numel(mn_files);
figure()

yyaxis('left')
plot(x .* [1;1], (mn_files + [1 -1] .* std_files./sqrt(n_obs_files))', 'r-', x, mn_files, '*-');
ylabel('Mean')

yyaxis('right')
plot(x, std_files)
ylabel('STD')

title(ttl)
xticks(x)
xticklabels(cellfun(@(x) x(4:5), f', 'uni', 0))
xlabel('Subset size')

end

function plot_detection_rates(summary, f)

summary.dr = [summary.dr; [1 0]];
f = [f; {'full '}];

mn = summary.dr(:, 1);
ci = mn + [1 -1] .* summary.dr(:, 2);

x = 1:numel(mn);

plot(x .* [1; 1], ci', 'r-', x, mn, '*')
xticklabels(cellfun(@(x) x(1:5), f', 'uni', 0))
title('Detection rates')
ylabel('Proportion');

end