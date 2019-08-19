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