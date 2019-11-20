function [wave_fit, mea] = wave_prop(mea, metric, varargin)
% Inputs: 
%	mea:		Matfile or struct with seizure epoch
%	dataToFit:  which type of data to fit the plane to (delays, maxdescent, events) 
% 
% Name-value parameter pairs:
%	fitMethod:  use 'bos' or 'nyc' method of fitting a plane to the data
%					(default: 'bos')
%	showPlots:  display detailed plots of data
%   T:          Time window to use for computing coherence in 'delays'
%					data. 
%					(default: 10) [s]
%   halfWin:	Time window surrounding firing discharges for computing
%					'maxdescent' and 'events'. 
%					(default: 50) [ms]

%% Parse inputs
p = inputParser;
allMetrics = {...
	'delays', 'maxdescent', 'events', ...
	'deviance', 'rising', 'falling'};
validMetrics = @(x) any(validatestring(x, allMetrics));

addRequired(p, 'mea', @(x) isstruct(x) || strcmpi(class(x), 'matlab.io.MatFile'));
addRequired(p, 'metric', validMetrics);

parse(p, mea, metric, varargin{:})
struct2var(p.Results)

if ~isfield(mea, 'params'), mea.params = init_mea_params(); end

fit_method = mea.params.fit_method;
half_win = mea.params.half_win;
exclude = mea.params.exclude;
thresh = mea.params.thresh;
T = mea.params.T;
band = mea.params.delay_band;
show_plots = mea.params.show_plots;

%% Convert mea to struct if it is not writable
if ~isstruct(mea)
	if ~exist(mea.Properties.Source, 'file')
		error('File not found');
	elseif ~mea.Properties.Writable
		mea = load(mea.Properties.Source);
		mea.Time = mea.Time();
	end
end

if exclude  % exclude non-spiking channels
	mea = exclude_channels(mea);
end

%% Assign wave fitting method
switch lower(fit_method)
	case 'bos'
		fit_wave = @(data, position) fit_wave_bos(data, position, metric);
	case 'nyc'
		fit_wave = @(data, position) ...
			fit_wave_nyc(data, position, metric);
end

%% Compute fits
[wave_fit, mea] = compute_waves(mea, fit_wave, show_plots, metric, ...
	T, half_win, thresh, band);

end

function [wave_fit, mea] = compute_waves(mea, fit_wave, show_plots, ...
    metric, T, half_win, thresh, band)
%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency


%% Pull variables from mea
position = mea.Position;

try
	position(mea.BadChannels, :) = [];
catch ME
	switch ME.message
		case 'Reference to non-existent field ''BadChannels''.'
			disp('No bad channels found.')
		otherwise
			rethrow(ME)
	end
end

Time = mea.Time;
Time = Time();
POS = (position - min(position)) ./ ... 
	arrayfun(@(ii) min(diff(unique(position(:, ii)))), [1 2]) + 1;% set to integers 1 .. n

% Load method specific variables
switch metric
	case {'rising'; 'falling'; 'deviance'}
		[computeTimes, mea] = get_waveTimes(mea);
		[lfp, skipfactor, mea] = get_lfp(mea);
% 		lfp = lfp - mean(lfp, 2);
		TimeMs = downsample(Time, skipfactor) * 1e3;
		lfp = lfp ./ std(lfp(TimeMs < 0, :));
        if strcmpi(metric, 'rising')
            dir = 1;
            if thresh == Inf, thresh = 10; end
        else
            dir = -1;
            if thresh == Inf, thresh = 20; end
        end
	case 'maxdescent'
		[computeTimes, mea] = get_waveTimes(mea);
		[lfp, skipfactor, mea] = get_lfp(mea);
		TimeMs = downsample(Time, skipfactor) * 1e3;
	case 'events'
		[computeTimes, mea] = get_waveTimes(mea);                          % Get discharge times
		if ~any(strcmpi(properties(mea), 'event_inds'))                    % If event times aren't already computed
			[~, ~, mea] = mua_events(mea);                                 % ... compute them
		end
		[nT, numCh] = size(mea.mua);
		[timeInds, chInds] = ind2sub([nT, numCh], mea.event_inds);         % Get the indices and channels of each spike
		TimeMs = Time * 1000;                                              % Convert times to ms
		spike_times = TimeMs(timeInds);                                    % Get the spike times in ms
        warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId
	case 'delays'
		[computeTimes, mea] = get_waveTimes(mea);  % Get discharge times
		
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
		Time = downsample(Time, skipfactor);
		mea.lfp = lfp; mea.skipfactor = skipfactor;
		TimeMs = Time * 1e3;
% 		compute_inds = arrayfun(@(t) find([TimeMs(:); Inf] >= t, 1), computeTimes);
		[params, ~] = set_coherence_params(Time, T, band);
		plotTitles = strrep(mea.Name, '_', ' ');
% 		computeTimes = TimeMs(compute_inds);
% 		samplingRate = mea.SamplingRate / skipfactor;
		[~, center] = min(sum((position - mean(position)).^2, 2));         % find the most central electrode
end
		
% Sizing variables
numCh = length(position);
numWaves = numel(computeTimes);

assignin('base', 'mea', mea);
%% Open a video file if PLOT is set to true

Name = sprintf('%s_wave_prop_%s', mea.Name, metric);
if show_plots
	v = VideoWriter(Name);
	plotTitles = [strrep(mea.Name, '_', ' ') ' (' metric ')'];
% 	v.FrameRate = 50;
	open(v); 
	h = figure; fullwidth(true);
	
	addy = sub2ind(max(POS), POS(:, 1), POS(:, 2));
end

%% Estimate wave direction at each discharge time

% Initialize arrays
beta = nan(3, numWaves);  % fit parameters
V = nan(2, numWaves);     % wave velocity (psuedo-inverse of beta)
p = nan(1, numWaves);     % certainty
[dataout, positionout] = deal(cell(numWaves, 1));  % wave passage times

for ii = 1:numWaves  % estimate wave velocity for each discharge
	position = POS;
	t = computeTimes(ii);
	if show_plots, img = nan(max(POS)); end
	switch metric
		case 'events'
			
% 			inds = logical((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
			inds = logical((spike_times >= t - half_win) .* ...            % Get spike times within halfWin ms
				(spike_times <= t + half_win));                            % ... of discharge at time t
% 			dataToPlot = nan(numCh, 1);
			pos_inds = chInds(inds);
            st = spike_times(inds);
            temp = arrayfun(@(ii) mean(st(pos_inds == ii)), 1:numCh);  % Find the mean spike time on each channel
			dataToPlot = temp - (t - half_win) + 1;  % imagesc mean spike time on each channel
% 			[~, pos_inds, data] = find(spike_times(inds, :));              % 
			data = spike_times(inds);
			temp = data(:)';
			position = position(pos_inds, :);
            if numel(data) == 0, continue, end
		
		case {'rising'; 'falling'; 'deviance'}
			inds = (TimeMs >= t - half_win) & (TimeMs <= t + half_win);  % Select the window around the event time
			temp = (smoothdata(lfp(inds, :), 'movmean', 5));  % A little smoothing to get rid of artefacts
			temp = (temp - temp(1, :));  % Set initial value as baseline
            if thresh == -Inf, threshI = quantile(max(dir*temp), .25) / 2; 
            else, threshI = thresh; end
			data = arrayfun(@(ii) ...  % Find where each channel deviates thresh from baseline
				find([dir * (temp(:, ii)); threshI] - threshI >= 0, 1), 1:size(temp, 2));
			data(data > size(temp, 1)) = size(temp, 1);
			data = data(:);
			tt = TimeMs(inds);
            dataToPlot = data;
			dataToPlot(~isnan(data)) = ...
                tt(data(~isnan(data))) - (t - half_win);  % Convert to time (ms)
			pos_inds = 1:numCh;
            if numel(unique(data)) < 3, continue, end
			
		case 'maxdescent'
			
			inds = logical((TimeMs >= t - half_win) .* (TimeMs <= t + half_win));  % Select the window around the event time
			temp = (smoothdata(lfp(inds, :), 'movmean', 5));  % A little smoothing to get rid of artefacts
			temp = temp - temp(1, :);  % set first time point as baseline (for visualization early)
			
			[change, time_point] = min(diff(temp, 1, 1));  % Find time of maximal descent
			non_decreasing = change >= 0;
            tt = TimeMs(inds);
			dataToPlot = tt(time_point(:)) - (t - half_win);
			dataToPlot(non_decreasing) = nan;
			
			pos_inds = 1:numCh;
			data = time_point;
			pos_inds(non_decreasing) = [];
			data(non_decreasing) = [];
			position(non_decreasing, :) = [];
			
			if numel(unique(data)) < 3, continue, end

		case 'delays'

			inds = and((TimeMs >= t - T / 2 * 1e3), ...
				(TimeMs <= t + T / 2 * 1e3));  % Select the window around the event time
			params.T = range(TimeMs(inds)) / 1e3;
% 			inds = compute_inds(ii) : (compute_inds(ii) + T * samplingRate - 1);
			if inds(end) > size(lfp, 1)
				break;
			end
			fprintf('Estimating waves at time %d/%d\n', ii, numWaves)
			temp = lfp(inds, :);
			[coh, phi, freq, coh_conf] = ...
				compute_coherence(temp, params, 'pairs', center);          % compute the coherence over the selected interval
			[delay, ~, ~] = compute_delay(coh, coh_conf, -phi, freq);      % compute delays on each electrode based on coherence
			data = 1e3 * delay(center,:)';                                 % we use delays relative to the center (converted to ms)
			dataToPlot = data + half_win;
			pos_inds = 1:numCh;
    end
    
    % Use the largest cluster of data points
    dataS = sort(data(:));
    diffdata = diff(dataS);  % calculate gaps between data points
    gaps = [0; ...
        find(diffdata > max(2*std(diffdata, 'omitnan'), 4)); ...
        numel(data)];  % find large gaps between datapoints
    [~, cM] = max(diff(gaps));  % find the largest group without a gap
    bounds = dataS(gaps(cM:cM+1) + [1; 0]);
    data(data < bounds(1) | data > bounds(2)) = nan;
    
    dataout{ii} = data;
    positionout{ii} = position;
	[beta(:, ii), V(:, ii), p(ii)] = fit_wave(data, position);
	
	if isnan(p(ii))
		% don't waste time (risk errors) plotting if a wave can't be fit
		continue
	end
	if show_plots
		figure(h); clf
		[p1, p2] = plot_wave_fit(position, data, beta(:, ii));
		title(p1, sprintf('%s\n %0.3f s', plotTitles, t / 1e3));
		title(p2, sprintf('p=%.2g', p(ii)))
		
% 		figure(h(2));
		img(addy) = dataToPlot;
		subplot(236); 
		p3 = imagesc(img', [-1 2*half_win]); axis xy
        colormap(h, 1 - parula(2 * half_win + 2));
		xlabel('X'); ylabel('Y');
		colorbar();
		cmap = h.Colormap;
		cInds = round(dataToPlot) + 1;
		cInds = min(max(cInds, 1), 2*half_win+2);
		if strcmpi(metric, 'events'), cInds = cInds(pos_inds); end
		
		subplot(2,3,4:5);
		plot_details(temp, pos_inds, cmap, cInds, metric); 
		axis tight; grid on;
		title(sprintf('%s\n %0.3f s', plotTitles, t / 1e3));
		if any(strcmpi(metric, {'maxdescent', 'deviance', 'rising', 'falling'}))
			valid = ~isnan(data);
			hold on; plot(dataToPlot(valid), ...
                temp(sub2ind(size(temp), ...
                data(valid)', pos_inds(valid))), 'r*'); hold off
		end
		try
            frame1 = getframe(h);
        catch ME
            disp(ii)
            disp(ME);
            disp(ME.stack);
            continue
        end
		writeVideo(v, frame1)
		
	end
	
end
Z = angle([1 1i] * V);
% Zu = unwrap(Z);

if show_plots, close(v); end

wave_fit.beta = beta;
wave_fit.V = V;
wave_fit.p = p;
wave_fit.Z = Z;
% wave_fit.Zu = Zu;
wave_fit.computeTimes = computeTimes;
wave_fit.Name = Name;
wave_fit.data = dataout;
wave_fit.position = positionout;


mea.wave_fit = wave_fit;

end

function [beta, V, p] = fit_wave_bos(data, position, metric)
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
	
	try  % If lfp is in mea, the load, otherwise filter
		lfp = mea.lfp;
		if any(strcmpi(fieldnames(mea), 'skipfactor'))
			skipfactor = mea.skipfactor;
		else
			skipfactor = 1;
		end

	catch ME
		if ~strcmp(ME.identifier, 'MATLAB:nonExistentField')
			rethrow(ME)
		end
		lfp = filter_mea(mea, {'lfp'});
		skipfactor = lfp.skipfactor;
		lfp = lfp.lfp;
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
