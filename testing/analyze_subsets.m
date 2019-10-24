data=load('testing/subsets/CUCX4_2_M.mat', '*stats');

f = fieldnames(data);
nF = numel(f);

% Convenience
mn =@(varargin) mean(varargin{:}, 'omitnan');
sd =@(varargin) std(varargin{:}, 'omitnan');
ci =@(data, outliers) sd(data) ./ sqrt(sum(~outliers));

[summary.dr, ...
	summary.dspeed, ...
	summary.dtheta] = deal(nan(nF, 2));

for ii = 1:nF
	struct2var(data.(f{dii}))
	
	% Detection rates
	outliers = isoutlier(detection_rate);
	detection_rate(outliers) = nan;
	summary.dr(ii, 1) = mn(detection_rate);
	summary.dr(ii, 2) = sd(detection_rate) ./ sqrt(sum(~outliers));
	
	% Differences in speed
	outliers = isoutlier(dspeed, 2);
	dspeed(outliers) = nan;
	dspeed_mn = mn(dspeed, 2);
	summary.dspeed(ii, 1) = mn(dspeed_mn);
	summary.dspeed(ii, 2) = ci(dspeed_mn, isoutlier(dspeed_mn));
	
	% Differences in direction
	
	
	
end