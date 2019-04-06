function [wave_fit, mea] = wave_prop(mea, method, PLOT, varargin)

if ~exist('PLOT', 'var')
	PLOT = false;
end
if ~isstruct(mea)
	if ~mea.Properties.Writable
		mea = load(mea.Properties.Source);
	end
end
	
switch lower(method)
	case 'bos'
		[wave_fit, mea] = bos_method(mea, PLOT, varargin);
	case 'nyc'
		[wave_fit, mea] = nyc_method(mea, PLOT);
end

end

function [wave_fit, mea] = nyc_method(mea, PLOT, varargin)
%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency

METHOD = 'lfp';

for i = 1:numel(varargin)
	switch varargin{i}
		case 'method'
			METHOD = varargin{i + 1};
		otherwise
			disp(varargin{i});
			error('Argument not recognized.');
	end
end

position = mea.Position;
try
	position(mea.BadChannels, :) = [];
catch ME
	disp(ME)
	disp('No bad channels found.')
end

if PLOT
	v = VideoWriter([mea.Name '_wave_prop_nyc']);
% 	v2 = VideoWriter([mea.Name '_wave_prop_nyc_2']);
	Name = strrep(mea.Name, '_', ' ');
% 	v.FrameRate = 30;
	open(v); 
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
	mea.lfp = lfp;
	mea.skipfactor = skipfactor;
end


%% Find discharge times

try 
	waveTimes = mea.waveTimes; 
catch ME
	disp(ME);
	disp('Computing wave times.');
	if ~isstruct(mea), mea = load(mea.Properties.Source); end
	[waveTimes, mea] = get_discharge_times(mea);
	mea.waveTimes = waveTimes;
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
numCh = length(position);
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
		frame1 = getframe(h(1));
		
		figure(h(2));
		img(addy) = [data{:}];
		subplot(122); imagesc(img); axis xy
		colorbar();
		cmap = parula(range([data{:}]) + 1);
		subplot(121); set(gca, 'ColorOrder', cmap(round([data{:}] - min([data{:}]) + 1), :), 'NextPlot', 'replacechildren'); 
		plot(temp); axis tight
		title(sprintf('%s\n %0.3f s', Name, t / 1e3));
		hold on; plot([data{:}]', temp(sub2ind(size(temp), [data{:}], 1:numCh))', '*'); hold off
		frame2 = getframe(h(2));
		
		frame.cdata = [frame1.cdata; frame2.cdata];
		frame.colormap = [];
		writeVideo(v, frame)
		
	end
	
end
Z = angle(complex(V(1, :), V(2, :)));
Zu = unwrap(Z);

if PLOT, close(v); end

wave_fit.beta = beta;
wave_fit.V = V;
wave_fit.p = p;
wave_fit.Z = Z;
wave_fit.Zu = Zu;
try
	mea.wave_fit_nyc = wave_fit;
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
src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);
psig = zeros(N, 1);
delays = zeros(N, length(position), length(position), 'single');
if PLOT
	phis = cell(N, 1);
end

for i = 1:N  % For each interval during the seizure
	
	% get time indices over which to compute coherence
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * ds_freq - 1);
	
	% compute the coherence over the selected interval
	fprintf('Estimating waves at time %d/%d\n', i, N)
	try
	[~, center] = min((position(:,1) - mean(position(:,1))).^2 +...         % find the most central electrode
				  (position(:,2) - mean(position(:,2))).^2);
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params, 'pairs', center);
	% compute delays on each electrode based on coherence
	[delay, ~, ~] = compute_delay(coh, coh_conf, phi, freq);
	
	
	delay = delay(center,:);                                                % we use delays relative to the center

	delays(i, :) = delay;
	
	if PLOT
		phi(coh < coh_conf) = NaN;
		phis{i} = phi;
	end
	
	% fit plane to delays
	[src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :), psig(i)] = ...
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
