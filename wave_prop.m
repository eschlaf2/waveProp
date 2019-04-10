function [wave_fit, mea] = wave_prop(mea, dataToFit, varargin)
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

allDataToFit = {'delays', 'maxdescent', 'events'};
allFitMethods = {'nyc', 'bos'};

validate = @(x, all) any(validatestring(x, all));

addRequired(p, 'mea', @(x) isstruct(x) || strcmpi(class(x), 'matlab.io.MatFile'));
addRequired(p, 'dataToFit', @(x) validate(x, allDataToFit));
addParameter(p, 'fitMethod', 'bos', @(x) validate(x, allFitMethods));
addParameter(p, 'showPlots', true, @islogical);
addParameter(p, 'T', 10, @isnumeric);
addParameter(p, 'halfWin', 50, @isnumeric);

parse(p, mea, dataToFit, varargin{:})
struct2var(p.Results)

%% Convert mea to struct if it is not writable
if ~isstruct(mea)
	if ~exist(mea.Properties.Source, 'file')
		error('File not found');
	elseif ~mea.Properties.Writable
		mea = load(mea.Properties.Source);
		mea.Time = mea.Time();
	end
end

%% Assign wave fitting method
switch lower(fitMethod)
	case 'bos'
		fit_wave = @(data, position) fit_wave_bos(data, position, dataToFit);
	case 'nyc'
		fit_wave = @(data, position) ...
			fit_wave_nyc(data, position, dataToFit);
end

%% Compute fits
[wave_fit, mea] = compute_waves(mea, fit_wave, showPlots, dataToFit, ...
	T, halfWin);

end

function [wave_fit, mea] = compute_waves(mea, fit_wave, showPlots, dataToFit, T, halfWin)
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
	disp(ME)
	disp('No bad channels found.')
end

Time = mea.Time;
POS = position;

% Load method specific variables
switch dataToFit
	case 'maxdescent'
		[computeTimes, mea] = get_waveTimes(mea);
		[lfp, skipfactor, mea] = get_lfp(mea);
		TimeMs = downsample(Time, skipfactor) * 1e3;
	case 'events'
		[computeTimes, mea] = get_waveTimes(mea);
		if ~any(strcmpi(properties(mea), 'event_inds'))
			[~, ~, mea] = mua_events(mea);
		end
		% Create a sparse array of spike times
		[nT, numCh] = size(mea.mua);
		spike_times = sparse(nT, numCh);
		spike_times(mea.event_inds) = 1;
		TimeMs = Time * 1000;  % Convert times to ms
		spike_times = TimeMs' .* spike_times;
	case 'delays'
		
		[lfp, skipfactor, mea] = get_lfp(mea);
		Time = downsample(Time, skipfactor);
		TimeMs = Time * 1e3;
		[params, compute_inds] = set_coherence_params(mea, Time, T);
		plotTitles = strrep(mea.Name, '_', ' ');
		computeTimes = TimeMs(compute_inds);
		samplingRate = mea.SamplingRate / skipfactor;
		[~, center] = min(sum((position - mean(position)).^2, 2));        % find the most central electrode
end
		
% Sizing variables
numCh = length(position);
numWaves = numel(computeTimes);

assignin('base', 'mea', mea);
%% Open a video file if PLOT is set to true

Name = sprintf('%s_wave_prop_%s', mea.Name, dataToFit);
if showPlots
	v = VideoWriter(Name);
	plotTitles = [strrep(mea.Name, '_', ' ') ' (' dataToFit ')'];
% 	v.FrameRate = 30;
	open(v); 
	h = figure; fullwidth(true);
	
	addy = sub2ind([10 10], position(:, 1), position(:, 2));
end

%% Estimate wave direction at each discharge time

% Initialize arrays
beta = nan(3, numWaves);  % fit parameters
V = nan(2, numWaves);     % wave velocity (psuedo-inverse of beta)
p = nan(1, numWaves);     % certainty

for ii = 1:numWaves  % estimate wave velocity for each discharge
	t = computeTimes(ii);
	if showPlots, img = nan(10); end
	switch dataToFit
		case 'events'
			
			inds = logical((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
			dataToPlot = nan(numCh, 1);
			temp = mean(spike_times(inds, :));
			[~, ch, temp] = find(temp);
% 			[~, so] = sort(ch);
			dataToPlot(ch) = temp;
			[~, pos_inds, data] = find(spike_times(inds, :));
			temp = data(:)';
			position = POS(pos_inds, :);
			
		case 'maxdescent'
			
			inds = logical((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
			temp = smoothdata(lfp(inds, :), 'movmedian', 10);                                           % Pull out window around wave
			[~, data] = min(diff(temp, 1, 1));                             % Find time of maximal descent
			data = data(:);
% 			temp = diff(temp, 1, 1);
			dataToPlot = data;
			pos_inds = 1:numCh;


		case 'delays'

			inds = compute_inds(ii) : (compute_inds(ii) + T * samplingRate - 1);
			fprintf('Estimating waves at time %d/%d\n', ii, numWaves)
			temp = lfp(inds, :);
			[coh, phi, freq, coh_conf] = ...
				compute_coherence(temp, params, 'pairs', center);  % compute the coherence over the selected interval
			[delay, ~, ~] = compute_delay(coh, coh_conf, phi, freq);       % compute delays on each electrode based on coherence
			data = delay(center,:)';                                       % we use delays relative to the center
			dataToPlot = data;
			pos_inds = 1:numCh;

	end
		
	[beta(:, ii), V(:, ii), p(ii)] = fit_wave(data, position);
	
	if showPlots
		figure(h);
		[p1, p2] = plot_wave_fit(position, data, beta(:, ii));
		title(p1, sprintf('%s\n %0.3f s', plotTitles, t / 1e3));
		title(p2, sprintf('p=%.2g', p(ii)))
		
% 		figure(h(2));
		img(addy) = dataToPlot;
		subplot(236); 
		p3 = imagesc(img); axis xy
		colorbar();
		cmap = h.Colormap;
		cInds = round((dataToPlot - min(dataToPlot))/range(dataToPlot) * (length(cmap) - 1)) + 1;
		if strcmpi(dataToFit, 'events'), cInds = cInds(pos_inds); end
		subplot(2,3,4:5);
		plot_details(temp, pos_inds, cmap, cInds, dataToFit); 
		axis tight; grid on;
		title(sprintf('%s\n %0.3f s', plotTitles, t / 1e3));
		if strcmpi(dataToFit, 'maxdescent')
			hold on; plot(dataToPlot, temp(sub2ind(size(temp), dataToPlot', 1:numCh)), 'r*'); hold off
		end
		frame1 = getframe(h);
% 		frame2 = getframe(h(2));
		
% 		frame.cdata = [frame1.cdata; frame2.cdata];
% 		frame.colormap = [];
		writeVideo(v, frame1)
		
	end
	
end
Z = angle(complex(V(1, :), V(2, :)));
Zu = unwrap(Z);

if showPlots, close(v); end

wave_fit.beta = beta;
wave_fit.V = V;
wave_fit.p = p;
wave_fit.Z = Z;
wave_fit.Zu = Zu;
wave_fit.computeTimes = computeTimes;
wave_fit.Name = Name;

try
	mea.wave_fit = wave_fit;
catch ME
	disp(ME)
	disp('Fit not saved.')
end

end


function [wave_fit, mea] = bos_method(mea, PLOT, varargin)


%% Parse input and set defaults

BAND = [1 13];                  % Select a frequency range to analyze
T = 10;                         % Length of recording (s)
W = 2;                          % Bandwidth
ntapers = 2*(T * W) - 1;        % Choose the # of tapers.
OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)
Name = strrep(mea.Name, '_', ' ');

if any(cellfun(@(v) strcmpi(v, 'T'), varargin{1}))
	ind = find(cellfun(@(v) strcmp(v, 'T'), varargin{1}));
	T = varargin{1}{ind + 1};
end

if PLOT
	% 
	PLOT = 'plot';
	figure;
	img_dir = sprintf('%s_wave_prop_bos_T%d', mea.Name, T);
	if ~exist(img_dir, 'dir')
		mkdir(img_dir);
	end
else
	PLOT = '';
end

try
	lfp = mea.lfp;  % import lfp band
	if any(strcmpi(fieldnames(mea), 'skipfactor'))
		skipfactor = mea.skipfactor;
	else
		skipfactor = 1;
	end
catch ME
	lfp = filter_mea(mea, [], 'lfp');
	skipfactor = lfp.skipfactor;
	lfp = lfp.lfp;
	
end

position = mea.Position;
if any(strcmpi(fieldnames(mea), 'BadChannels'))
	position(mea.BadChannels, :) = [];
end

[~, center] = min((position(:,1) - mean(position(:,1))).^2 +...         % find the most central electrode
			  (position(:,2) - mean(position(:,2))).^2);

%% Set parameters
% Note: the chronux toolbox must be on the path

% Downsample lfp
	
ds_freq = 1e3;  % downsampled frequency (Hz)
samplingRate = mea.SamplingRate;
lfp_freq = samplingRate / skipfactor;
ds_step = floor(max(1, lfp_freq / ds_freq));  % how to step through the data
lfp = lfp(1 : ds_step : end, :);  % downsample
		
ds_freq = lfp_freq / ds_step;  % recompute DS_FREQ in case of rounding
nsamp = size(lfp, 1);

time = mea.Time;  % import time
time = time();
time = time(1:skipfactor:end);  % downsample to match lfp
time = time(1 : ds_step : end);  % downsample further if necessary
% wave_times = mea.waveTimes;
padding = mea.Padding;

params.tapers = [T * W, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = ds_freq;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.


COMPUTE_INDS = [round(1 : OVERLAP_COMPLEMENT * ds_freq : nsamp - T * ds_freq) nsamp];
% [~, COMPUTE_INDS] = arrayfun(@(t) min(abs(t - time)), wave_times / 1e3);

TS = find(time(COMPUTE_INDS) > 0, 1);
TE = find(time(COMPUTE_INDS) > (time(end) - max(padding(2), T)), 1) - 1;
if isempty(TE)
	TE = numel(COMPUTE_INDS);
end

COMPUTE_INDS = COMPUTE_INDS(TS : TE);  % only compute while in seizure
N = numel(COMPUTE_INDS);


%%

% Initialize arrays
beta = zeros(N, 3);
src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);
psig = zeros(N, 1);
delays = zeros(N, length(position), 'single');
if PLOT
	phis = cell(N, 1);
end

for i = 1:N  % For each interval during the seizure
	
	% get time indices over which to compute coherence
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * ds_freq - 1);
	
	% compute the coherence over the selected interval
	fprintf('Estimating waves at time %d/%d\n', i, N)
	try
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params, 'pairs', center);
	% compute delays on each electrode based on coherence
	[delay, ~, ~] = compute_delay(coh, coh_conf, phi, freq);
	
	delay = delay(center,:);                                                % we use delays relative to the center

	delays(i, :) = delay;
	
	if PLOT
		phi(coh < coh_conf) = NaN;
		phis{i} = squeeze(phi(center, :, :));
	end
	
	% fit plane to delays
	[beta(i, :), src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :), psig(i)] = ...
		estimate_wave(delay, position, PLOT);
	
	catch ME
		disp(ME)
		disp(i)
		continue
	end
	
	if PLOT
		title(sprintf('%s\n %0.3f s', Name, time(COMPUTE_INDS(i))));
		print(gcf, fullfile(img_dir, num2str(i, '%03d')), '-dpng')
% 		frame = getframe(gcf);
% 		writeVideo(v, frame);
	end

end
V = speed .* exp(1i * src_dir);  % wave velocity (for consistency with NYC method)
V = [real(V) imag(V)];

wave_fit = struct('Z', src_dir, 'V', V, 'sp', speed, ...
	'ci_Z', ci_dir, 'ci_sp', ci_sp,...
	'delays', delays, ...
	'freq', freq, ...
	'p', psig, ...							  % significance level
	'params', params, ...                     % analysis parameters
	'wave_times', time(COMPUTE_INDS) * 1e3);  % wave times in ms.
if PLOT
	wave_fit.phis = phis;
end
try
	mea.(sprintf('wave_fit_bos_T%d', T)) = wave_fit;
catch ME
	disp(ME);
	disp('Fit not saved to matfile!')
	[~, fn, ~] = fileparts(mea.Properties.Source);
	save([fn '_wave_fit'], wave_fit, '-v7.3');
end

if PLOT
	png2avi(img_dir);
end
end

function [beta, V, p] = fit_wave_bos(data, position, dataToFit)
	[beta, ~, ~, ~, ~, p] = estimate_wave(data, position);
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
		disp(ME);
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
		disp(ME);
		lfp = filter_mea(mea, [], {'lfp'});
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

function [params, compute_inds] = set_coherence_params(mea, Time, T)

	BAND = [1 13];                  % Select a frequency range to analyze
	W = 2;                          % Bandwidth
	NTAPERS = 2*(T * W) - 1;        % Choose the # of tapers.
	OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)
	
	samplingRate = mea.SamplingRate / mea.skipfactor;
	
	params.tapers = [T * W, NTAPERS];  % ... time-bandwidth product and tapers.
	params.Fs = samplingRate; % ... sampling rate
	params.pad = -1;                % ... no zero padding.
	params.fpass = BAND;            % ... freq range to pass
	params.err = [1 0.05];          % ... theoretical error bars, p=0.05.
	params.T = T;
	
% 	padding = mea.Padding;
	nsamp = numel(Time);                                % Total number of samples
	step = OVERLAP_COMPLEMENT * samplingRate;           % Compute coherence every OVERLAP_COMPLEMENT
	lastSamplePoint = nsamp - T * samplingRate;         % Leave a long enough window at the end to calculate coherence
	compute_inds = round(1 : step : lastSamplePoint);    

% 	ts = find(Time(compute_inds) > 0, 1);
% 	te = find(Time(compute_inds) > (Time(end) - max(padding(2), T)), 1) - 1;
% 	if isempty(te)
% 		te = numel(compute_inds);
% 	end
% 
% 	compute_inds = compute_inds(ts : te);  % only compute while in seizure
		
end

function [] = plot_details(data, position, cmap, cInds, dataToFit)
	
	switch dataToFit
		case 'events'
			[cInds, so] = sort(cInds);
			set(gca, 'ColorOrder', cmap(cInds, :), 'NextPlot', 'replacechildren'); 
			for ii = 1:numel(data)
				plot(data(so(ii)), position(ii), '*'); hold on;
			end
		otherwise
			nanmask = isnan(cInds);
			cInds(nanmask) = [];
			data = data(:, ~nanmask);
			set(gca, 'ColorOrder', cmap(cInds, :), 'NextPlot', 'replacechildren'); 
			plot(data)
	end
end

function [p1, p2] = plot_wave_fit(position, data, beta)

	X = position(:, 1);
	Y = position(:, 2);
	beta = beta(:);
	

	p1 = subplot(2,3,1:2);
	scatter3(p1, X, Y, data, 200, data, 'filled');hold on;        

	[XX, YY] = meshgrid(sort(unique(X)),sort(unique(Y)));
	Z = double(position * beta(1:2) + beta(end));
	f = scatteredInterpolant(X, Y, Z);
	ZZ = f(XX,YY);
	mesh(p1, XX, YY, ZZ, ...
		'FaceColor', 'interp', 'FaceAlpha', 0.8, 'LineStyle', 'none') %interpolated

	legend(p1, 'Data','Regression fit');
	xlabel('cm');ylabel('cm');zlabel('Second');
	
	hold off;
	% Plot the projection along the velocity axis
	p2 = subplot(233);
	P_v_axis = position * beta(1:2) / norm(beta(1:2));
	plot(p2, P_v_axis, data, '.');
	hold on;
	plot(p2, P_v_axis, Z);
	title('Projection along the velocity vector');
	xlabel('cm');
	ylabel('Second'); colormap(p2, 'hot')
	hold off;
	

end 
