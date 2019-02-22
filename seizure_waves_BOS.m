function [mea] = seizure_waves_BOS(mea, PLOT)

%% Parse input
if ~exist('PLOT', 'var')
	if ~feature('ShowFigureWindows')
		PLOT = true;
	else
		PLOT = false;
	end
end

% data = mea.lfp;
position = [mea.X mea.Y];

%% Set parameters
% Note: the chronux toolbox must be on the path

% Downsample lfp
DS_FREQ = 1e3;  % downsampled frequency (Hz)
DS_STEP = floor(mea.SamplingRate / DS_FREQ);  % how to step through the data
DS_FREQ = mea.SamplingRate / DS_STEP;  % recompute DS_FREQ in case of rounding

lfp = mea.lfp;  % import lfp band
lfp = lfp(1 : DS_STEP : end, :);  % downsample
NSAMP = size(lfp, 1);

time = mea.Time;  % import time
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

src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);

for i = 1:N
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * DS_FREQ - 1);
	
	fprintf('Estimating waves at time %d/%d\n', i, N)
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params);
	
	[delay, delay_ci_lo, delay_ci_up] = compute_delay(coh, coh_conf, phi, freq);
	
	[src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :)] = estimate_wave(delay, position, '');

end

BOS = struct('src_dir', src_dir, 'speed', speed, 'ci_dir', ci_dir, 'ci_sp', ci_sp);
mea.BOS = BOS;

