function wave_fit = wave_prop(mea, method)


switch lower(method)
	case 'bos'
		wave_fit = bos_method(mea);
	case 'nyc'
		wave_fit = nyc_method(mea);
end

end

function wave_fit = nyc_method(mea)
%% Compute wave propagation at each discharge time as described in 
% Liou, Jyun You, et al. ?Multivariate Regression Methods for Estimating
% Velocity of Ictal Discharges from Human Microelectrode Recordings.?
% Journal of Neural Engineering, vol. 14, no. 4, NIH Public Access, 2017,
% p. 044001, doi:10.1088/1741-2552/aa68a6. 
% 
% All time units should be converted to ms for consistency

halfWin = 50;  % half window around discharge event (ms)
Time = mea.Time;
Time = Time();

%% Find discharge times

fr = mea.firingRate;
mask = mean(fr) >= 1/60;  % exclude channels with mean firing rate less than one spike per minute (Liou et al., 2018) ?\cite{Liou2018a}
meanFr = mean(fr(:, mask), 2);

[~, waveTimes] = findpeaks(meanFr, ...  % find peaks in mean firing rate
	mea.SamplingRate / 1e3, ...  % ... in ms 
	'minpeakprom', 100 * std(diff(meanFr)), ...  % ... use discrete peaks
	'minpeakdistance', 100);  % ... peaks should be at least 100 ms apart

padding = mea.Padding;
waveTimes = waveTimes - padding(1) * 1e3;  % Account for padding (in ms)

mea.waveTimes = waveTimes;


% Create an array of spike times
T = nan(size(mea.mua), 'single');
T(mea.event_inds) = 1;
TimeMs = Time * 1000;  % Convert times to ms
T = TimeMs' .* T;
% T(T == 0) = nan;

% Sizing variables
numCh = numel(mea.X);
numWaves = numel(waveTimes);


%% Estimate wave direction at each discharge time

% Initialize arrays
beta = nan(3, numWaves);  % fit parameters
V = nan(2, numWaves);  % wave velocity (psuedo-inverse of beta)
p = nan(1, numWaves);  % certainty
position = [mea.X, mea.Y];

for i = 1:numWaves  % estimate wave velocity for each discharge
	t = waveTimes(i);
	inds = find((TimeMs >= t - halfWin) .* (TimeMs <= t + halfWin));
	events = T(inds, :)';
	events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
	for ch = 1:numCh
		temp = events{ch};
		temp(isnan(temp)) = [];
		events{ch} = temp;
	end
	[beta(:, i), V(:, i), p(i)] = ...
		SpatialLinearRegression(events, position, ...
		'switch_plot', 0, 'Lossfun','L2');
end
Z = angle(complex(V(1, :), V(2, :)));
Zu = unwrap(Z);

wave_fit.beta = beta;
wave_fit.V = V;
wave_fit.p = p;
wave_fit.Z = Z;
wave_fit.Zu = Zu;
mea.wave_fit_nyc = wave_fit;

end


function wave_fit = bos_method(mea)

position = [mea.X mea.Y];

%% Set parameters
% Note: the chronux toolbox must be on the path

% Downsample lfp
DS_FREQ = 1e3;  % downsampled frequency (Hz)
DS_STEP = floor(mea.SamplingRate / DS_FREQ);  % how to step through the data
DS_FREQ = mea.SamplingRate / DS_STEP;  % recompute DS_FREQ in case of rounding

try
	lfp = mea.lfp;  % import lfp band
catch ME
	lfp = filter_mea(mea, [], 'lfp');
	lfp = lfp.lfp;
end
lfp = lfp(1 : DS_STEP : end, :);  % downsample
NSAMP = size(lfp, 1);

time = mea.Time;  % import time
time = time();
padding = mea.Padding;
time = time(1 : DS_STEP : end);  % downsample


BAND = [1 13];                  % Select a frequency range to analyze
T = 2;							% Length of recording (s)
W = 2;                          % Bandwidth
ntapers = 2*(T * W) - 1;        % Choose the # of tapers.
params.tapers = [T * W, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = DS_FREQ;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.
OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)

COMPUTE_INDS = [round(1 : OVERLAP_COMPLEMENT * DS_FREQ : NSAMP - T * DS_FREQ) NSAMP];
TS = find(time(COMPUTE_INDS) > 0, 1);
TE = find(time(COMPUTE_INDS) > (time(end) - padding(2)), 1) - 1;

COMPUTE_INDS = COMPUTE_INDS(TS : TE - T);  % only compute while in seizure
N = numel(COMPUTE_INDS);
%%

% Initialize arrays
src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);
psig = zeros(N, 1);

for i = 1:N  % For each interval during the seizure
	
	% get time indices over which to compute coherence
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * DS_FREQ - 1);
	
	% compute the coherence over the selected interval
	fprintf('Estimating waves at time %d/%d\n', i, N)
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params);
	
	% compute delays on each electrode based on coherence
	[delay, delay_ci_lo, delay_ci_up] = compute_delay(coh, coh_conf, phi, freq);
	
	% fit plane to delays
	[src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :), psig(i)] = ...
		estimate_wave(delay, position, '');

end
V = speed .* exp(1i * src_dir);  % wave velocity (for consistency with NYC method)
V = [real(V) imag(V)];

wave_fit = struct('Z', src_dir, 'V', V, 'sp', speed, ...
	'ci_Z', ci_dir, 'ci_sp', ci_sp,...
	'p', psig, ...							  % significance level
	'params', params, ...                     % analysis parameters
	'wave_times', time(COMPUTE_INDS) * 1e3);  % wave times in ms.
mea.wave_fit_bos = wave_fit;

end