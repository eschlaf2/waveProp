function wave_fit = wave_prop(mea, method, PLOT, varargin)

if ~exist('PLOT', 'var')
	PLOT = false;
end
switch lower(method)
	case 'bos'
		wave_fit = bos_method(mea, PLOT, varargin);
	case 'nyc'
		wave_fit = nyc_method(mea, PLOT);
end

end

function wave_fit = nyc_method(mea, PLOT)
%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency

METHOD = 'lfp';
% DOWNSAMPLE = true;

position = [mea.X, mea.Y];

if PLOT
	v1 = VideoWriter([mea.Name '_wave_prop_nyc_1']);
	v2 = VideoWriter([mea.Name '_wave_prop_nyc_2']);
	Name = strrep(mea.Name, '_', ' ');
% 	v.FrameRate = 30;
	open(v1); open (v2);
	h(1) = figure; fullwidth();
	h(2) = figure; fullwidth();
	
	img = nan(10);
	addy = sub2ind([10 10], position(:, 1), position(:, 2));
end

halfWin = 50;  % half window around discharge event (ms)
Time = mea.Time;
Time = Time();
try
	lfp = mea.lfp;
catch ME
	lfp = filter_mea(mea, [], {'lfp'});
	lfp = lfp.lfp;
end
if ~isprop(mea, 'skipfactor')
	mea.skipfactor = 1;
end
skipfactor = mea.skipfactor;

%% Find discharge times

try 
	waveTimes = mea.waveTimes; 
catch ME
	waveTimes = get_discharge_times(mea);
end

if strcmpi(METHOD, 'events')
	% Create an array of spike times
	T = nan(size(mea.mua), 'single');
	T(mea.event_inds) = 1;
	TimeMs = Time * 1000;  % Convert times to ms
	T = TimeMs' .* T;
else
	TimeMs = downsample(Time, skipfactor) * 1e3;
end

% Sizing variables
numCh = numel(mea.X);
numWaves = numel(waveTimes);


%% Estimate wave direction at each discharge time

% Initialize arrays
beta = nan(3, numWaves);  % fit parameters
V = nan(2, numWaves);  % wave velocity (psuedo-inverse of beta)
p = nan(1, numWaves);  % certainty

for i = 1:numWaves  % estimate wave velocity for each discharge
	t = waveTimes(i);
	inds = find((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
	switch METHOD
		case 'events'
			events = T(inds, :)';
			events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
			for ch = 1:numCh
				temp = events{ch};
				temp(isnan(temp)) = [];
				events{ch} = temp;
			end
			data = events;
		case 'lfp'
			temp = lfp(inds, :);
			[~, data] = min(diff(temp, 1, 1));  % maximal descent
			data = num2cell(data);
	end
	
	if PLOT, figure(h(1)); end
	[beta(:, i), V(:, i), p(i)] = ...
		SpatialLinearRegression(data, position, ...
		'Lossfun', 'L2', 'switch_plot', PLOT);
	if PLOT
		title(sprintf('%s\n %0.3f s', Name, t / 1e3));
		frame = getframe(h(1));
		writeVideo(v1, frame);
		
		figure(h(2));
		img(addy) = [data{:}];
		subplot(122); imagesc(img); axis xy
		colorbar();
		cmap = parula(range([data{:}]) + 1);
		subplot(121); set(gca, 'ColorOrder', cmap(round([data{:}] - min([data{:}]) + 1), :), 'NextPlot', 'replacechildren'); 
		plot(temp); axis tight
		title(sprintf('%s\n %0.3f s', Name, t / 1e3));
		hold on; plot([data{:}]', temp(sub2ind(size(temp), [data{:}], 1:numCh))', '*'); hold off
		frame = getframe(h(2));
		writeVideo(v2, frame)
		
	end
	
end
Z = angle(complex(V(1, :), V(2, :)));
Zu = unwrap(Z);

if PLOT, close(v1); close(v2); end

wave_fit.beta = beta;
wave_fit.V = V;
wave_fit.p = p;
wave_fit.Z = Z;
wave_fit.Zu = Zu;
mea.wave_fit_nyc = wave_fit;

end


function wave_fit = bos_method(mea, PLOT, varargin)


%% Parse input and set defaults

BAND = [1 13];                  % Select a frequency range to analyze
T = 10;		% Length of recording (s)
W = 2;                          % Bandwidth
ntapers = 2*(T * W) - 1;        % Choose the # of tapers.
OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)

if any(cellfun(@(v) strcmpi(v, 'T'), varargin{1}))
	ind = find(cellfun(@(v) strcmp(v, 'T'), varargin{1}));
	T = varargin{1}{ind + 1};
end

if PLOT
	PLOT = 'plot';
	figure;
	v = VideoWriter([mea.Name '_wave_prop_bos']);
	Name = strrep(mea.Name, '_', ' ');
% 	v.FrameRate = 30;
	open(v);
else
	PLOT = '';
end

position = [mea.X mea.Y];

%% Set parameters
% Note: the chronux toolbox must be on the path

try
	lfp = mea.lfp;  % import lfp band
catch ME
	lfp = filter_mea(mea, [], 'lfp');
	lfp = lfp.lfp;
end

% Downsample lfp
if isprop(mea, 'skipfactor')
	DS_STEP = mea.skipfactor;
else
	DS_FREQ = 1e3;  % downsampled frequency (Hz)
	DS_STEP = floor(mea.SamplingRate / DS_FREQ);  % how to step through the data
	lfp = lfp(1 : DS_STEP : end, :);  % downsample

end
DS_FREQ = mea.SamplingRate / DS_STEP;  % recompute DS_FREQ in case of rounding
NSAMP = size(lfp, 1);

time = mea.Time;  % import time
time = time();
time = time(1 : DS_STEP : end);  % downsample
% wave_times = mea.waveTimes;
padding = mea.Padding;

params.tapers = [T * W, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = DS_FREQ;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.


COMPUTE_INDS = [round(1 : OVERLAP_COMPLEMENT * DS_FREQ : NSAMP - T * DS_FREQ) NSAMP];
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
src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);
psig = zeros(N, 1);

parfor i = 1:N  % For each interval during the seizure
	
	% get time indices over which to compute coherence
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * DS_FREQ - 1);
	
	% compute the coherence over the selected interval
	fprintf('Estimating waves at time %d/%d\n', i, N)
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params);
	
	% compute delays on each electrode based on coherence
	[delay, ~, ~] = compute_delay(coh, coh_conf, phi, freq);
	
	% fit plane to delays
	[src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :), psig(i)] = ...
		estimate_wave(delay, position, PLOT);
	if PLOT
		title(sprintf('%s\n %0.3f s', Name, time(COMPUTE_INDS(i))));
		frame = getframe(gcf);
		writeVideo(v, frame);
	end

end
V = speed .* exp(1i * src_dir);  % wave velocity (for consistency with NYC method)
V = [real(V) imag(V)];

wave_fit = struct('Z', src_dir, 'V', V, 'sp', speed, ...
	'ci_Z', ci_dir, 'ci_sp', ci_sp,...
	'p', psig, ...							  % significance level
	'params', params, ...                     % analysis parameters
	'wave_times', time(COMPUTE_INDS) * 1e3);  % wave times in ms.
try
	mea.wave_fit_bos = wave_fit;
catch ME
	warning(ME);
end

if PLOT
	close(v);
end
end
