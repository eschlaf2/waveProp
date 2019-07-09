% delays_by_localregression
if isinteger(phi); phi = single(phi) / 1e4; end
fband = [25 50];
finds = (f >= fband(1)) & (f <= fband(2));
phi(C <= confC) = nan;
dphi = permute(phi, [2 1 3]);  % permute array to f x t x pairs
dphi = gradient(dphi, f(2) - f(1));  % d/df

winsz = 3;  % Hz
dim = 1; % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
degree = 2;  % linear
method = 'rlowess';

delays = -matlab.internal.math.localRegression(dphi, winsz, dim, ...
            nanflag, degree, method, f);
        
%% Imagesc delays

clims = quantile(delays(:), [.025 .975]);
for ii = 1:10
    temp = delays(:, :, ii);
    imagesc(t, f, temp, clims);
    axis xy; colorbar; ylim(fband); drawnow(); pause();
end 

%% Fit delays

predictors = [ones(size(f)); f; f.^2; f.^3]';
predictors = [ones(size(f))]';
delaysR = reshape(delays, length(f), []);
[polyfit, ~, ~, ~, stats] = arrayfun(@(ii) regress(delaysR(finds, ii), predictors(finds, :)), 1:length(delaysR), 'uni', 0);
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
metrics = {'delays', 'delays_all', 'delays_T02_fband25_50'};
res = plot_wave_polar(res, metrics, sig, ax(1), ax(2));
legend(metrics)




