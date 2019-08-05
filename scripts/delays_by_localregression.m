% delays_by_localregression

winsz = 3;  % Hz
thresh = 5e-2; 
df = f(2) - f(1); 
fband = [1 50];
% MASK = false;

if isinteger(phi); phi = -single(phi) / 1e4; end
finds = (f >= fband(1)) & (f <= fband(2));

transform = @(A) reshape(permute(A(:, finds, :), [2 1 3]), sum(finds), []);
Cf = transform(C);  % limit to band of interest
phif = transform(phi);  % ... same for phi
% phif = smoothdata(unwrap(phif), 1, 'rlowess', winsz / df);  % ... unwrap
% dphi = padarray(diff(unwrap(phif)), 1, 'pre');
[~, dphi] = gradient(unwrap(phif), df);  % Compute the gradient of phi wrt freq
dphi(Cf <= confC) = nan;  % set insignificant values to nan

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

	dphi(~mask) = nan;
end

% dphi = diff(phif) / df;

%% Delays
dim = 1;  % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
degree = 1;  % constant
method = 'movmed';

% delays = matlab.internal.math.localRegression(dphi, winsz / df, dim, ...
%             nanflag, degree, method, f);
[delays, wn] = smoothdata(dphi, dim, method, winsz / df, nanflag);	

delaysR = reshape(delays, size(delays, 1), numel(t), []);

%% Imagesc delays

if 0
    
clims = quantile(delays(:), [.025 .975]);
for ii = 1:10
	inds = (ii - 1) * numel(t) + 1: ii * numel(t);
    temp = delays(:, inds);
    imagesc(t, f(finds), temp, clims);
	line(t, ones(size(t)) * 13, 'color', 'green', 'linewidth', 2);
    axis xy; colorbar; ylim(fband); drawnow(); pause();
end 

end
%% Imagesc delaysR

if 0
    
for ii = 1:10
    imagesc(t, f(finds), delaysR(:, :, ii), clims);
    line(t, ones(size(t)) * 13, 'color', 'green', 'linewidth', 2);
    axis xy; colorbar; ylim(fband); drawnow(); pause();
end

end

%% All wave directions

MIN_RATIO_FINITE = .25;
[nf, nt, np] = size(delaysR);
[Z, pdel, pct] = deal(nan(nf, nt));
pos = position(pairs(:, 2), :);

warning('off', 'stats:statrobustfit:IterationLimit');
for ii = 1:nf
    for jj = 1:nt
        delays2fit = squeeze(delaysR(ii, jj, :));
        finite = sum(isfinite(delays2fit));
		pct(ii, jj) = finite;
        if finite <= max(MIN_RATIO_FINITE * size(pos, 1), 3); continue; end  % check enough delay data is not NaN.
        [beta,stats] = robustfit(pos, delays2fit, 'fair');                     % fit the delay vs two-dimensional positions
        H = [0 1 0; 0 0 1];  
        c = [0; 0];
        ptemp = linhyptest(beta, stats.covb, c, H, stats.dfe);  
        if ptemp > thresh || isnan(ptemp); continue; end
        V = pinv(beta(2:3));
        Z(ii, jj) = angle([1 1i] * V(:));
		pdel(ii, jj) = ptemp;
    end
end

%% Imagesc Z (angles computed using delays)
[tt, ff] = ndgrid(t, f(finds));
cmap = hsv(80);
h = pcolor(tt, ff, Z'); h.LineStyle = 'none'; colormap(cmap); colorbar;
h.Parent.CLim = [-pi pi];
line(t, 13 * ones(size(t)), 'color', 'black', 'linewidth', 2)
xlabel('Time (s)');
ylabel('Freq (Hz)')
title(sprintf('%s Seizure %d', pat, seizure));
% ylim([fband(1) fband(2)-winsz]);

%% Compare to measure
% compareto = 'events';

if 0
	
wavefit = load(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop_1.mat', pat, seizure), compareto);
comparemetric = interp1(wavefit.(compareto).computeTimes, wavefit.(compareto).Z, t * 1e3, 'nearest');
diffs = (comparemetric - Z)';
diffs(diffs > pi) = diffs(diffs > pi) - 2 * pi;
diffs(diffs < -pi) = diffs(diffs < -pi) + 2 * pi;
[tt, ff] = ndgrid(t, f(finds));
cmap = [flipud(hot(40)); hot(40)];
h = pcolor(tt, ff, diffs); h.LineStyle = 'none'; colormap(cmap); colorbar;
h.Parent.CLim = [-pi pi] * 1.1;
h.Parent.Color = .5*[1 1 1];
line(t, 13 * ones(size(t)), 'color', 'black', 'linewidth', 2)
xlabel('Time (s)');
ylabel('Freq (Hz)')
title(sprintf('%s Seizure %d\n%s', pat, seizure, compareto));

end

%% Fit delays

if 0

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

end

%% Fit waves

if 0
    
delays2fit = bfit.o0;
% delaystofit = -bfit.o1;
[beta, ~, ~, ~, ~, pdel] = arrayfun(@(ii)...
    estimate_wave(delays2fit(ii, :), position(pairs(:, 2), :)), ...
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

end
%% for comparison with other delay algorithms

if 0
    
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

end


%% print plots

if 0
    
files = dir('*_Neuroport_10_10_cohgram_ds_T02.mat');
% files = files(29:end);
counter = 1;

close all
figure(1); fullwidth(true);
for f = files'
    fname = f.name
    info = strsplit(fname, '_');         
    pat = info{1}; seizure = str2double(info{2}(8:end));
    load(fname)
    subplot(2, 2, mod(counter - 1, 4) + 1);
    delays_by_localregression;
    if ~mod(counter, 4)
        print(1, sprintf('delays_by_lr_%s_%02d', compareto, counter), '-dpng')
        clf;
    end
    counter = counter + 1;
    disp(counter)
    
end
print(1, sprintf('delays_by_lr_%d', counter), '-dpng')
clf;

end