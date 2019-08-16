% delays_by_localregression
% basename = compute_coherograms(pat, seizure, T, W, DS, units, toi);  %
% (See quickscript)
load(basename, ...  % Load the following from basename
    'f', 'params', 'phi', 'C', 't', 'confC', 'pat', ...
    'seizure', 'position', 'pairs')

winsz = 3;  % Hz
thresh = 5e-2; 
df = diff(f(1:2)); 
fband = params.fpass;
if ~exist('toi', 'var'); toi = [-Inf Inf]; end
if ~exist('type', 'var'); type = 'group'; end
% MASK = false;

if isinteger(phi); phi = -single(phi) / 1e4; end
finds = (f >= fband(1)) & (f <= fband(2));
tinds = (t >= toi(1)) & (t <= toi(2));
t = t(tinds); f = f(finds);

[nt, nf, np] = size(C);

transform = @(A) reshape(permute(A(tinds, finds, :), [2 1 3]), sum(finds), []);
Cf = transform(C);  % limit to band of interest
phif = transform(phi);  % ... same for phi
clear C phi
% phif = smoothdata(unwrap(phif), 1, 'rlowess', winsz / df);  % ... unwrap
% dphi = padarray(diff(unwrap(phif)), 1, 'pre');

disp('Computing gradient')
[~, dphi] = gradient(unwrap(phif), df);  % Compute the gradient of phi wrt freq

% dphi = padarray(diff(unwrap(phif)), 1, 'pre') / df;
dphi(Cf <= confC) = nan;  % set insignificant values to nan

if exist('MASK', 'var') && MASK  % Keep only the longest streak of significant data points
	dphi = keep_streak(dphi, winsz / df);
end


%% Delays
disp('Computing delays')

dim = 1;  % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
degree = 1;  % constant
method = 'movmed';

% delays = matlab.internal.math.localRegression(dphi, winsz / df, dim, ...
%             nanflag, degree, method, f);
switch type
	case 'group'
		delays = smoothdata(dphi, dim, method, winsz / df, nanflag);
	case 'phase'
		delays = smoothdata(phif ./ f', dim, method, winsz / df, nanflag);
end

delaysR = reshape(delays, size(delays, 1), numel(t), []);
clear delays phif

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
disp('Computing directions')

MIN_RATIO_FINITE = .25;
pos = position(pairs(:, 2), :);

% p = gcp();
warning('off', 'stats:statrobustfit:IterationLimit');
H = [0 1 0; 0 0 1];  
c = [0; 0];

packages = ver;
arrayfun_is_faster = ...  % Use parfor if possible
    any(cellfun(@(n) strcmpi(n, 'parallel computing toolbox'), {packages.Name}));

if ~arrayfun_is_faster
    disp('Using arrayfun ...')
    tic %#ok<*UNRCH>
    delaysR2 = reshape(delaysR, [], np)';  % reshape the delays so that each column is a single time-freq point and rows are pairs
    num_unique = arrayfun(@(ii) ...  % Number of unique, non-nan values in each column
        numel(unique(delaysR2(isfinite(delaysR2(:, ii)), ii))), 1:length(delaysR2));
    finite = sum(isfinite(delaysR2));  % number of non-nan, non-inf delays in each column
    Z = nan * num_unique;  % Initialize Z with nans
    inds = ...  % Need at least 3 unique values and MIN_RATIO_FINITE finite values
        find((finite >= max(MIN_RATIO_FINITE * np, 3)) & (num_unique >= 3));  
    
    disp('Computing fits')
    tic
    [beta, stats] = ...  % Compute fits
        arrayfun(@(ii) robustfit(pos, delaysR2(:, ii), 'fair'), ...
        inds, 'uni', 0);  
    toc
    
    disp('Computing significance')
    tic
    pdel = ...  % ... and significance
        cellfun(@(A, B) linhyptest(A, B.covb, c, H, B.dfe), beta, stats);  
    toc
    
    indsB = find(pdel < thresh);  % Keep significant fits
    disp('Computing Z')
    tic
    Z(inds(indsB)) = ...  % Fill in wave directions for good fits
        arrayfun(@(ii) angle(pinv(beta{ii}([2 3])) * [1; 1i]), indsB);
    Z = reshape(Z, nf, nt);
    toc
    
    toc
else
    disp('Using parfor ...')
    tic
    [Z, pdel, pct] = deal(nan(nf, nt));
    parfor ii = 1:nf  % For each frequency
        if ~mod(ii, 100), fprintf('ii=%d/%d\n', ii, nf), end
        for jj = 1:nt  % ... and time point
            delays2fit = squeeze(delaysR(ii, jj, :));  % ... collect the delays for each pair
            num_unique = numel(unique(delays2fit(isfinite(delays2fit))));  % ... calculate the number of finite unique values
            if num_unique < 3, continue; end  % check that there are enough unique points to fit a plane
            finite = sum(isfinite(delays2fit));  % Calculate the number of finite delays
            pct(ii, jj) = finite / np;  % (FYI only)
            if finite <= max(MIN_RATIO_FINITE * np, 3); continue; end  % check enough delay data is not NaN.
            [beta,stats] = robustfit(pos, delays2fit, 'fair');  % fit the delay vs two-dimensional positions
            ptemp = linhyptest(beta, stats.covb, c, H, stats.dfe);  % Compute significance of fit
            if ptemp >= thresh || isnan(ptemp); continue; end  % Stop if fit is not significant
            V = pinv(beta(2:3));  % velocity is the pseudoinvervse of the fitted wave
            Z(ii, jj) = angle([1 1i] * V(:));  % Z is the angle of the velocity
            pdel(ii, jj) = ptemp;  % Store p-value
        end
    end
    toc
end

%% Imagesc Z (angles computed using delays)
emilys_pcolor(t, f * units, Z_gpdelays', 'cmap', hsv(80), 'clim', [-pi,pi]);
line(t, 13 * ones(size(t)), 'color', 'black', 'linewidth', 2)
xlabel('Time (s)');
ylabel('Freq (Hz)')
title(sprintf('%s Seizure %d\%s delays', pat, seizure, type));

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
	predictors = ones(size(f))';
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