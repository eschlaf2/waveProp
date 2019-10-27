data=load('testing/subsets/CUCX4_2_M.mat', '*stats');

f = fieldnames(data);
nF = numel(f);
dim = 1;  % dim=1: returns stats for each wave time; 
          % dim=2: returns states for each trial

% Convenience
mn =@(data) mean(data, dim, 'omitnan');
sd =@(data) std(data, [], dim, 'omitnan');
ci =@(data) 2 * sd(data, [], dim) ./ sqrt(sum(isfinite(data), dim));

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
	dspeed_mn = mn(dspeed, dim);
	summary.dspeed(ii, 1, :) = dspeed_mn;
	summary.dspeed(ii, 2, :) = ci(dspeed_mn);
	
	% Differences in direction
	outliers = isoutlier(abs(dtheta), dim);
	dtheta(outliers) = nan;
	summary.dtheta(ii, 1, :) = circ_mean(dtheta, [], dim, 'omitnan');
	confmean = circ_confmean(dtheta, [], dim, 'omitnan');
	if ~isreal(confmean), confmean = angle(confmean); end
	summary.dtheta(ii, 2, :) = confmean;
	
	
end