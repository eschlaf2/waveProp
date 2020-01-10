function [wave_fit, mea] = wave_prop(mea, metric, varargin)
% Inputs: 
%	mea:		Matfile or struct with seizure epoch
%	metric:  which type of data to fit the plane to (delays, maxdescent, events) 
% 
% Name-value parameter pairs:
%	fitMethod:  use 'bos' or 'nyc' method of fitting a plane to the data
%					(default: 'bos')
%	showPlots:  display detailed plots of data
%   T:          Time window to use for computing coherence in 'delays'
%					data. 
%					(default: 10) [s]
%   half_win:	Time window surrounding firing discharges for computing
%					'maxdescent' and 'events'. 
%					(default: 50) [ms]

%% Parse inputs

% Main
[mea, metric] = parse_inputs(mea, metric, varargin{:});

% Parameters
fit_method = mea.params.fit_method;
half_win = mea.params.half_win;
exclude = mea.params.exclude;
thresh = mea.params.thresh;
T = mea.params.T;
band = mea.params.delay_band;
show_plots = mea.params.show_plots;

%%
% exclude non-spiking channels
if exclude, mea = exclude_channels(mea); end

%% Assign wave fitting method
switch lower(fit_method)
	case 'bos'
		fit_wave = @(data, position) fit_wave_bos(data, position, metric);
	case 'nyc'
		fit_wave = @(data, position) ...
			fit_wave_nyc(data, position, metric);
end

%% Compute fits
% [wave_fit, mea] = compute_waves(mea, fit_wave, show_plots, metric, ...
% 	T, half_win, thresh, band);
% 
% 	
% end
% 
% function [wave_fit, mea] = compute_waves(mea, fit_wave, show_plots, ...
%     metric, T, half_win, thresh, band)
%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency


%% Set static loop variables

NAME = sprintf('%s_wave_prop_%s', mea.Name, metric);
[COMPUTE_TIMES, mea] = get_waveTimes(mea);
[TIME_MS, POS, LFP, RF, E, D, PLT] = set_globs;

% Convenience variables
NUM_WAVES = numel(COMPUTE_TIMES);

assignin('base', 'mea', mea);

%% Estimate wave direction at each discharge time

% Initialize arrays
beta = nan(3, NUM_WAVES);  % fit parameters
V = nan(2, NUM_WAVES);     % wave velocity (psuedo-inverse of beta)
p = nan(1, NUM_WAVES);     % certainty
[data_out, position_out] = deal(cell(NUM_WAVES, 1));  % wave passage times

for ii = 1:NUM_WAVES  % estimate wave velocity for each discharge
	
	t = COMPUTE_TIMES(ii);
	inds = ...
		get_time_inds;
	pos_inds = ...
		which_positions;
	[fit_data, time_series_data, im_data] = ...
		get_data;
	
	if numel(unique(fit_data)) < 3, continue, end
	
	reduced_data = ...
		use_largest_cluster;
	
    data_out{ii} = reduced_data;
    position_out{ii} = POS(pos_inds, :);
	
	[beta(:, ii), V(:, ii), p(ii)] = ...
		fit_wave(reduced_data, POS(pos_inds, :));
	
	% don't waste time (risk errors) plotting if a wave can't be fit
% 	if isnan(p(ii)), continue, end
	
	generate_plots;
	
end

compile_results;

	
%% Nested functions

% Setup
function [time_ms, position, lfp, RF, E, D, plt] = set_globs
	
	% Common variables
	position = mea.Position;
	if isfield(mea, 'BadChannels'), position(mea.BadChannels, :) = []; end
	time = mea.Time();
	time_ms = time * 1e3;
	% POS = (position - min(position)) ./ ... 
	% 	arrayfun(@(ii) min(diff(unique(position(:, ii)))), [1 2]) + 1;  % set to integers 1 .. n

	% Method specific variables
	lfp = [];  % maxdescent, rising, falling, deviance, delays
	RF.thresh = thresh;  % rising, falling
	E = [];  % events variables
	D = [];  % delays variables
	plt = [];  % plotting variables
	zscore =@(X) (X - nanmean(X)) ./ nanstd(X);
	if ismember(lower(metric), {'maxdescent', 'rising', 'falling', 'deviance'})
		[lfp, skipfactor, mea] = get_lfp(mea);
		lfp = zscore(smoothdata(lfp, 'movmean', 5));  % A little smoothing to get rid of artefacts
		time_ms = downsample(time_ms, skipfactor);
	end
	
	% Method specific variables
	switch metric
		case {'rising'; 'falling'; 'deviance'}
			lfp = lfp ./ std(lfp(time_ms < 0, :));
			if strcmpi(metric, 'rising')
				RF.dir = 1;
				if thresh == Inf, RF.thresh = 10; end
			else
				RF.dir = -1;
				if thresh == Inf, RF.thresh = 20; end
			end
		case 'events'
			if ~isfield(mea, 'event_inds')  % If event times aren't already computed
				[~, ~, mea] = mua_events(mea);  % ... compute them
			end
			if ~isfield(mea, 'event_mat_size'), mea.event_mat_size = size(mea.mua); end
			[time_idx, E.ch_idx] = ind2sub(mea.event_mat_size, mea.event_inds);         % Get the indices and channels of each spike
			E.spike_times = time_ms(time_idx);                                    % Get the spike times in ms
			warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
			
		case 'delays'
			% Get band appropriate lfp
			filterband = ceil(band + [-.5 .5] .* band);
			if filterband(2) >= round(mea.SamplingRate / 2)
				filterband(2) = round(mea.SamplingRate / 2) - 1; 
			end
			b = fir1(150, 2 * filterband(2) / mea.SamplingRate);  % lo-pass to just over upper band
			lfp = single(filtfilt(b, 1, double(mea.Data)));
			skipfactor = round(mea.SamplingRate / max(1e3, 2*filterband(2)));
			lfp = downsample(lfp, skipfactor);
			lfp(:, mea.BadChannels) = [];
			[mea.lfp, mea.skipfactor] = deal(lfp, skipfactor);
			
			time_ms = downsample(time_ms, skipfactor);
			
			[D.params, ~] = set_coherence_params(time_ms / 1e3, T, band);
			[~, D.center] = min(sum((position - mean(position)).^2, 2));         % find the most central electrode
	end
	
	% Plotting variables
	if show_plots
		plt.v = VideoWriter(NAME);
		plt.plot_titles = [strrep(mea.Name, '_', ' ') ' (' metric ')'];
		open(plt.v); 
		plt.h = figure; fullwidth(true);
		plt.addy = sub2ind(max(position), position(:, 1), position(:, 2));
	end
	
end


% Analysis
function inds = get_time_inds

	switch metric
		case 'events'
			% Get spike times within half_in ms
			inds = (E.spike_times >= t - half_win) & (E.spike_times <= t + half_win);  % ... of discharge at time t                            

		case {'rising'; 'falling'; 'deviance'}
			inds = (TIME_MS >= t - half_win) & (TIME_MS <= t + half_win);  % Select the window around the event time
		case 'maxdescent'
			% Select time points within half_win ms of the event time
			inds = logical((TIME_MS >= t - half_win) & (TIME_MS <= t + half_win));  

		case 'delays'
			% Select time_points within T/2 s of the event time
			inds = (TIME_MS >= t - T / 2 * 1e3) & (TIME_MS <= t + T / 2 * 1e3);  
	end
end

function pos_inds = which_positions

	if strcmpi(metric, 'events')
		pos_inds = E.ch_idx(inds);
	else
		pos_inds = 1:length(POS);
	end

end

function [fit_data, time_series_data, im_data] = get_data
% Returns data to be fit depending on the metric (method)
% EVENTS
%     fit_data: spike times on each channel (though really this is
%         treated identically to the mean spike time on each channel.
%     im_data: mean spike time on each channel relative to window start
%     time_series_data: flattened fit_data to be plotted as a
%         raster
% RISING/FALLING/DEVIANCE
%     fit_data: time of threshold crossing in the rising, falling,
%         or agnostic (deviance) direction
%     im_data: same as fit_data
%     time_series_data: baseline adjusted lfp traces
% MAXDESCENT
%     fit_data: time of maximal descent in lfp (traces with max
%         descent on the boundary or traces that are non-decreasing
%         are excluded.
%     im_data: same as fit_data
%     time_series_data: baseline adjusted lfp traces
% DELAYS
%     fit_data: time delay (ms) of each electrode relative to the
%         central electrode
%     im_data: same as fit_data
%     time_series_data: lfp traces

	switch metric
		case 'events'
			
            fit_data = E.spike_times(inds);
			im_data = nan(length(POS), 1);
			[G, ID] = findgroups(pos_inds);
			im_data(ID) = splitapply(@mean, fit_data(:), G) - (t - half_win) + 1;
			time_series_data = reshape(fit_data(:), 1, []);
		
		case {'rising'; 'falling'; 'deviance'}
			
			time_series_data = LFP(inds, :);
			time_series_data = (time_series_data - time_series_data(1, :));  % Set initial value as baseline
			if RF.thresh == -Inf
				threshI = quantile(max(RF.dir * time_series_data), .25) / 2; 
			else
				threshI = thresh; 
			end
			tt = TIME_MS(inds);
			fit_data = arrayfun(@(ii) ...  % Find where each channel deviates thresh from baseline
				find([RF.dir * (time_series_data(:, ii)); threshI] - threshI >= 0, 1), 1:size(time_series_data, 2));
			fit_data(fit_data > size(time_series_data, 1)) = size(time_series_data, 1);
			fit_data = fit_data(:);
			im_data = tt(fit_data) - (t - half_win);    % Convert to time (ms) from window start
			
		case 'maxdescent'
			
			time_series_data = LFP(inds, :);
			time_series_data = time_series_data - time_series_data(1, :);  % set first time point as baseline (for visualization early)
			
			[change, time_point] = min(diff(time_series_data, 1, 1));  % Find time of maximal descent
			non_decreasing = change >= 0;  % Find non-decreasing traces
			bdry = (time_point == 1) | (time_point == size(time_series_data, 1) - 1);  % ... and traces with max descent on the boundary (these are often not part of the wave and confuse the analysis)
			inactive = range(time_series_data) < 2;
            tt = TIME_MS(inds);
			im_data = tt(time_point(:)) - (t - half_win);  % convert to time from window start
			im_data(non_decreasing | bdry | inactive) = nan;
			
			fit_data = time_point;
			fit_data(non_decreasing | bdry | inactive) = nan;
			
		case 'delays'

			D.params.T = range(TIME_MS(inds)) / 1e3;
			fprintf('Estimating waves at time %d/%d\n', ii, NUM_WAVES)
			time_series_data = LFP(inds, :);
			[coh, phi, freq, coh_conf] = ...
				compute_coherence(time_series_data, D.params, 'pairs', D.center);          % compute the coherence over the selected interval
			[delay, ~, ~] = compute_delay(coh, coh_conf, -phi, freq);      % compute delays on each electrode based on coherence
			fit_data = 1e3 * delay(D.center, :)';                                 % we use delays relative to the center (converted to ms)
			im_data = fit_data + half_win;
	end
end

function reduced_data = use_largest_cluster
    % Use the largest cluster of data points for lfp methods
% 	if ~ismember(metric, {'maxdescent', 'rising', 'falling', 'deviance'})
% 		reduced_data = fit_data; 
% 		return
% 	end

	dataS = sort(fit_data(:));
	if all(isnan(dataS)), reduced_data = fit_data; return; end
	dataS(isnan(dataS)) = [];  % excluded nan values
	diff_sorted = diff(dataS);  % calculate gaps between nearby data points
	big_gaps = isoutlier(unique(diff_sorted));  % find large gaps between datapoints
	divides = [0; find(big_gaps); numel(dataS)];  % divisions between clusters in sorted data
	cluster_sizes = diff(divides);  % size of clusters
	[~, largest] = max(cluster_sizes);  % choose the largest cluster

	bounds = dataS(divides(largest+[0 1]) + [1; 0]);
	reduced_data = fit_data;
	reduced_data(fit_data < bounds(1) | fit_data > bounds(2)) = nan;
end


% Output
function generate_plots
	if ~show_plots, return, end
	figure(PLT.h); clf
	img = nan(max(POS));
	
	% Top plots (wave fit and projection)
	[p1, p2] = plot_wave_fit(POS(pos_inds, :), reduced_data, beta(:, ii));
	title(p1, sprintf('%s\n %0.3f s', PLT.plot_titles, t / 1e3));
	title(p2, sprintf('p=%.2g', p(ii)))
	
	% Lower right plot (2D image)
	img(PLT.addy) = im_data;
	subplot(236); 
	imagesc(img', [-1 2*half_win]); axis xy
	colormap(PLT.h, 1 - parula(2 * half_win + 2));
	xlabel('X'); ylabel('Y');
	colorbar();
	cmap = PLT.h.Colormap;
	cInds = round(im_data) + 1;
	cInds = min(max(cInds, 1), 2*half_win+2);
	if strcmpi(metric, 'events'), cInds = cInds(pos_inds); end
	
	% Lower left plot (time series)
	subplot(2,3,4:5);
	plot_details(time_series_data, pos_inds, cmap, cInds, metric); 
	axis tight; grid on;
	title(sprintf('%s\n %0.3f s', PLT.plot_titles, t / 1e3));
	if any(strcmpi(metric, {'maxdescent', 'deviance', 'rising', 'falling'}))
		valid = ~isnan(reduced_data);
		hold on; plot(im_data(valid), ...
			time_series_data(sub2ind( ...
				size(time_series_data), ...
				fit_data(valid(:)), pos_inds(valid(:))) ...
			), 'r*'); 
		hold off
	end
	try
		frame1 = getframe(PLT.h);
	catch ME
		disp(ii)
		disp(ME);
		disp(ME.stack);
		return
	end
	writeVideo(PLT.v, frame1)

end

function compile_results
	
	if show_plots, close(PLT.v); end

	Z = angle([1 1i] * V);
	wave_fit.beta = beta;
	wave_fit.V = V;
	wave_fit.p = p;
	wave_fit.Z = Z;
	% wave_fit.Zu = Zu;
	wave_fit.computeTimes = COMPUTE_TIMES;
	wave_fit.Name = NAME;
	wave_fit.data = data_out;
	wave_fit.position = position_out;


	mea.wave_fit = wave_fit;
end

end

%% Helpers

function [beta, V, p] = fit_wave_bos(data, position, ~)
	[beta, ~, ~, ~, ~, p] = estimate_wave(data, position);
    
    % beta is invalid if nan or if the slope is 0 in both directions
    invalid = any(isnan(beta)) || all(beta(2:3).^2 < eps^2);
	if invalid
		V = nan;
		return
	end
	beta = circshift(beta, -1);
	V = pinv(beta(1:2));
end

function [] = fit_wave_nyc(data, position, dataToFit, showPlots)

% Convert to cell for SLR
	switch dataToFit
		case 'maxdescent'
			data = num2cell(data);
		case 'events'
			numCh = size(data, 2);
			data = mat2cell(data', ones(size(data, 1), 1), numCh);
			for ch = 1:numCh
				temp = data{ch};
				temp(isnan(temp)) = [];
				data{ch} = temp;
			end
		case 'delays'
			data = num2cell(data');
	end
	SpatialLinearRegression(data, position, ...
			'Lossfun', 'L2', 'switch_plot', showPlots);
end

function [waveTimes, mea] = get_waveTimes(mea)
	try % Find discharge times
		waveTimes = mea.waveTimes; 
	catch ME
		if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
			rethrow(ME);
		end
		disp('Computing wave times.');
		[waveTimes, mea] = get_discharge_times(mea);
	end
end

function [lfp, skipfactor, mea] = get_lfp(mea)
	
	new_samplingRate = 1e3;  % downsample to 1kHz
	
	if isfield(mea, 'lfp') && isfield(mea, 'skipfactor')
		[lfp, skipfactor] = deal(mea.lfp, mea.skipfactor);
	else
		lfp = filter_mea(mea, {'lfp'});
		[lfp, skipfactor] = deal(lfp.lfp, lfp.skipfactor);		
	end
	
	% Downsample to new_samplingRate
	samplingRate = mea.SamplingRate / skipfactor;
	new_skip = floor(samplingRate / new_samplingRate);
	if new_skip > 1
		lfp = downsample(lfp, new_skip);
		skipfactor = skipfactor * new_skip;
	end
	mea.lfp = lfp;
	mea.skipfactor = skipfactor;

end

function [params, compute_inds] = set_coherence_params(Time, T, band)

% 	band = [1 13];                  % Select a frequency range to analyze
	W = 3/T;                          % Bandwidth
	NTAPERS = 2*(T * W) - 1;        % Choose the # of tapers.
	OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)
	
	samplingRate = 1 / mean(diff(Time));
	
	params.tapers = [T * W, NTAPERS];  % ... time-bandwidth product and tapers.
	params.Fs = samplingRate; % ... sampling rate
	params.pad = 4;                 % ... 2^(ceil(log2(T)) + pad)
	params.fpass = band;            % ... freq range to pass
	params.err = [1 0.05];          % ... theoretical error bars, p=0.05.
	params.T = T;
	
	nsamp = numel(Time);                                % Total number of samples
	step = OVERLAP_COMPLEMENT * samplingRate;           % Compute coherence every OVERLAP_COMPLEMENT
	lastSamplePoint = nsamp - ceil(T * samplingRate) + 1;     % Leave a long enough window at the end to calculate coherence
	compute_inds = round(1 : step : lastSamplePoint);    
		
end

function [] = plot_details(data, position, cmap, cInds, metric)
	
	switch metric
		case 'events'
			[cInds, so] = sort(cInds);
			set(gca, 'ColorOrder', cmap(cInds, :), ...
                'NextPlot', 'replacechildren'); 
			for ii = 1:numel(data)
				plot(data(so(ii)), position(so(ii)), '*'); hold on;
			end
		otherwise
			nanmask = isnan(cInds);
			cInds(nanmask) = [];
			set(gca, 'ColorOrder', cmap(cInds, :), 'NextPlot', 'replacechildren'); 
			plot([1 size(data, 1)], [0 0], 'r'); hold on;
			data = data(:, ~nanmask);
			plot(data); hold off;
	end
end

function [p1, p2] = plot_wave_fit(position, data, beta)
    
	p1 = subplot(2,3,1:2);
    p2 = subplot(233);
    
	X = position(:, 1);
	Y = position(:, 2);
	beta = beta(:);
	

	scatter3(p1, X, Y, data, 200, data, 'filled'); hold(p1, 'on');        

    [XX, YY] = meshgrid(sort(unique(X)),sort(unique(Y)));
    Z = double(position * beta(1:2) + beta(end));
    if all(size(XX) > 2)
        f = scatteredInterpolant(X, Y, Z);
        ZZ = f(XX,YY);
        mesh(p1, XX, YY, ZZ, ...
            'FaceColor', 'interp', 'FaceAlpha', 0.8, ...
            'LineStyle', 'none') %interpolated

        legend(p1, 'Data','Regression fit');
    else
        legend(p1, 'Data')
    end
	xlabel('X (electrode)');ylabel('Y (electrode)');zlabel('Time (ms)');
	
	hold off;
	
    % Plot the projection along the velocity axis
	P_v_axis = position * beta(1:2) / norm(beta(1:2));
	plot(p2, P_v_axis, data, '.');
	hold on;
	plot(p2, P_v_axis, Z);
	title('Projection along the velocity vector');
	xlabel('cm');
	ylabel('Time (ms)'); colormap(p2, 'hot')
	hold off;
	

end 

function [mea, metric] = parse_inputs(mea, metric, varargin)
	p = inputParser;
allMetrics = {...
	'delays', 'maxdescent', 'events', ...
	'deviance', 'rising', 'falling'};
validMetrics = @(x) any(validatestring(x, allMetrics));

addRequired(p, 'mea', @(x) isstruct(x) || strcmpi(class(x), 'matlab.io.MatFile'));
addRequired(p, 'metric', validMetrics);

parse(p, mea, metric, varargin{:})
mea = p.Results.mea;
metric = p.Results.metric;

if ~isfield(mea, 'params'), mea.params = init_mea_params(); end

end
