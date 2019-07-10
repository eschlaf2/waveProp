% delays_by_localregression

winsz = 5;  % Hz
thresh = 5e-2; 
df = f(2) - f(1); 
fband = [1 50];
% MASK = false;

if isinteger(phi); phi = single(phi) / 1e4; end
finds = (f >= fband(1)) & (f <= fband(2));

transform = @(A) reshape(permute(A(:, finds, :), [2 1 3]), sum(finds), []);
Cf = transform(C);  % limit to band of interest
phif = unwrap(transform(phi));  % ... same for phi
phif(Cf <= confC) = nan;  % set insignificant values to nan

if exist('MASK', 'var') && MASK
	temp = padarray(Cf, [1 0], 'both');  % Pad with 0s
	starts = arrayfun(@(ii) ...  % find where streaks of significant coherence start
		find((temp(:, ii) <= confC) & (circshift(temp(:, ii), [-1 0]) > confC)), ...
		1:size(temp, 2), 'uni', 0);
	ends = arrayfun(@(ii) ...  % ... and where they end
		find((temp(:, ii) > confC) & (circshift(temp(:, ii), [-1 0]) <= confC)), ...
		1:size(temp, 2), 'uni', 0);
	[~, streak] = cellfun(@(s, e) max(e - s), starts, ends, 'uni', 0);  % Isolate the longest streak
	startinds = cellfun(@(s, ii) s(ii), starts, streak, 'uni', 0);  % Get the start ind for each streak
	endinds = cellfun(@(e, ii) e(ii), ends, streak, 'uni', 0);   % ... and the end ind
	mask = false(size(temp));  % Initialize a mask of longest streak of significant values
	for ii = 1:size(temp, 2)
		if endinds{ii} - startinds{ii} < winsz / df; continue; end
		mask(startinds{ii}:endinds{ii}, ii) = true;
	end
	mask([1 end], :) = [];  % unpad

	phif(~mask) = nan;
end

[~, dphi] = gradient(phif, df);  % Compute the gradient of phi wrt freq

%% Delays
dim = 1;  % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
% degree = 2;  % linear
method = 'rlowess';

% delays = -matlab.internal.math.localRegression(dphi, winsz, dim, ...
%             nanflag, degree, method, f);
[delays, wn] = smoothdata(-dphi, 1, method, winsz / df, nanflag);	
        
%% Imagesc delays

clims = quantile(delays(:), [.025 .975]);
for ii = 1:10
	inds = (ii - 1) * numel(t) + 1: ii * numel(t);
    temp = delays(:, inds);
    imagesc(t, f(finds), temp, clims);
	line(t, ones(size(t)) * 13, 'color', 'green', 'linewidth', 2);
    axis xy; colorbar; ylim(fband); drawnow(); pause();
end 

%% Fit delays

predictors = [ones(size(f)); f; f.^2; f.^3]';
predictors = [ones(size(f))]';
% delaysR = reshape(delays, length(f), []);
% delaysR = delaysR .* mask;
[polyfit, ~, ~, ~, stats] = arrayfun(@(ii) regress(delays(:, ii), predictors(finds, :)), 1:size(delays, 2), 'uni', 0);
polyfit = reshape(cat(2, polyfit{:}), size(predictors, 2), length(t), []);
polyfit = permute(polyfit, [2 3 1]);
pdel = reshape(cellfun(@(x) x(3), stats), length(t), []);

for ii = 1:size(polyfit, 3)
    fn = sprintf('o%d', ii - 1);
    bfit.(fn) = polyfit(:, :, ii); bfit.(fn)(pdel >= .05) = nan;
end

temp = bfit.o0;
clims = quantile(temp(:), [.01, .99]);
imagesc(t, pairs(:, 2), temp', clims); colorbar

% bfit.o0 = squeeze(median(delays, 1, 'omitnan'));

%% Fit waves
delaystofit = bfit.o0;
% delaystofit = -bfit.o1;
[beta, ~, ~, ~, ~, pdel] = arrayfun(@(ii)...
    estimate_wave(delaystofit(ii, :), position(pairs(:, 2), :)), ...
    1:numel(t), 'uni', 0);
pdel = [pdel{:}]';

% beta is invalid if nan or if the slope is 0 in both directions
invalid = cellfun(@(x) any(isnan(x)) || all(x(2:3).^2 < eps^2), beta);
beta = cellfun(@(x) circshift(x, -1), beta, 'uni', 0);
V = cell(length(t), 1);
V(~invalid) = cellfun(@(x) pinv(x(1:2)), beta(~invalid), 'uni', false);
V(invalid) = {[nan nan]};
V = cat(1, V{:});
beta(invalid) = {[nan; nan; nan]}; beta = [beta{:}]';
Z = angle(V * [1; 1i]);

%% for comparison with other delay algorithms
pat = 'c7';
plotnum = 1;
wave_dir_polarplots;
res.data.delays_all.Z = interp1(t, Z, res.time, 'nearest');
res.data.delays_all.p = interp1(t, pdel, res.time, 'nearest');
res.data.delays_all.beta = beta(interp1(t, 1:length(t), res.time, 'nearest', 'extrap'), :);
res.Z = [res.Z res.data.delays_all.Z];
res.p = [res.p res.data.delays_all.p];
metrics = {'events', 'delays_all', 'delays_T02_fband25_50'};
res = plot_wave_polar(res, metrics, sig, ax(1), ax(2));
legend(metrics)




