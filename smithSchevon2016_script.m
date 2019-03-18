% disp('Loading data')
% m = matfile('MG49_Seizure43.mat');
filename = m.Name;
outfile = [];
% outfile = [filename '_meaFilt'];

% The isfield() function does not return as expected for matfile objects so
% use this instead
isfield_or_var = @(struct_or_matfile, fieldname) ...
	any(strcmp(fieldnames(struct_or_matfile), fieldname));

%%

try
	mea = matfile([m.Name '_Neuroport_Filt']);
% 	ecog = m.ecog;
catch ME
	switch ME.identifier
		case 'MATLAB:MatFile:VariableNotInFile'
% 			errorLog = {errorLog; ME};
			clear ME
			if isempty(outfile)
				mea = m.Neuroport;
% 				ecog = m.ECoG;
			else
				mea = matfile(outfile, 'writable', true);
			end
		otherwise
			rethrow(ME)
	end
end
	
%% Make figures pres ready
defaultToPresentationFigs;
c = lines(7);
aic = @(predictors, dev) dev + 2 * size(predictors, 2);

tmax = Inf;
%% Filter

if ~isfield_or_var(mea, 'lfp')
	disp('Filtering');
% 	[mea, ecog] = filter_mea_ecog(mea, ecog);
	mea = filter_mea_ecog(mea, [], [], {'mua'});
	mea = filter_mea_ecog(mea, [], [], {'lfp'});
end

%% Get MUA events
% Negative peaks in this signal were detected and those peaks that exceeded 
% four times the s.d. of the signal in the negative direction, and that occurred 
% >1ms after immediately preceding peaks, were retained as multiunit timestamps 
% and waveforms (Smith et al., 2016)

disp('Getting MUA events');

if ~isfield_or_var(mea, 'event_inds')
	mea = mua_events(mea);
end
%% Compute firing rate
% Firing rate was calculated from these multiunit timestamps in 100-ms windows 
% every 25 ms (Smith et al., 2016)

disp('Computing firing rate');

if ~isfield_or_var(mea, 'firingRate')
	mea = mua_firing_rate(mea);
end
%% Test for recruitment
% Use the methods outlined in Schevon et al., 2012 to see if the MEA territory 
% was recruited into the seizure. This includes looking at the coherence and Fano 
% factor.

disp('Computing coherence')

if ~isfield_or_var(mea, 'wpli')
	mea = test_for_recruitment(mea, 'wpli');
end

disp('Computing Fano factor')

if ~isfield_or_var(mea, 'ff')
	mea = test_for_recruitment(mea, 'fano');
end
%% Compute wave velocity at discharges

disp('Performing spatial linear regression')

mea = wave_prop(mea);

%% Save figures

disp('Saving figures')

nn = strrep(m.Name, 'Seizure', '');

print(2, nn, '-dpng')
print(3, [nn '_FF'], '-dpng')
print(4, [nn, '_coh'], '-dpng')

%% Compute wave direction using BOS method

% Downsample lfp
position = [mea.X mea.Y];
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

PLOT = ''; % 'plot';
src_dir = zeros(N, 1);
speed = zeros(N, 1);
ci_dir = zeros(N, 2);
ci_sp = zeros(N, 2);

for i = 1:N
	inds = COMPUTE_INDS(i) : (COMPUTE_INDS(i) + T * DS_FREQ - 1);
	
	fprintf('Estimating waves at time %d/%d\n', i, N)
	[coh, phi, freq, coh_conf] = compute_coherence(lfp(inds, :), params);
	
	[delay, delay_ci_lo, delay_ci_up] = compute_delay(coh, coh_conf, phi, freq);
	
	[src_dir(i), speed(i), ci_dir(i, :), ci_sp(i, :)] = estimate_wave(delay, position, PLOT);

end

BOS = struct('src_dir', src_dir, 'speed', speed, 'ci_dir', ci_dir, 'ci_sp', ci_sp);
mea.BOS = BOS;

%% Define epochs and IDIs
% The ictal wavefront epoch was defined as the time of the first channel?s mean 
% minus its s.d. until the last channel?s mean plus its s.d. (Smith et al., 2016). 
% [I think it should read "max" so I used max instead of mean]

% disp('Defining epochs')
% firingRateSm = smoothdata(mea.firingRate, 'gaussian', 100);
% peakRate = zeros(size(mea.firingRate, 2), 1);
% frT = mea.Time;
% inds = logical((frT >= 0) .* (frT <= 20));  % Look in the first 20 seconds of the seizure
% frT = frT(inds);
% mea.epochs = zeros(3, 1);
% 
% for i = 1:numel(peakRate)
% 	[~, ind] = max((firingRateSm(inds, i)));
% 	peakRate(i) = frT(ind);
% end
% 
% peakRate([68 72]) = [];  % These two electrodes look very different... one is definitely some kind of inhibitory cell
% [firstT, firstCh] = min(peakRate);
% [lastT, lastCh] = max(peakRate);
% 
% firstT = firstT - std(peakRate);
% lastT = lastT + std(peakRate);
% 
% mea.epochs(1) = firstT;  % Recruitment start time
% mea.epochs(2) = lastT;  % Post recruitment start time
% 
% % Discharge times: find peaks in firing rate
% hga = (mean(mea.hga, 2) - mean(mean(mea.hga(mea.Time < 0, :), 2))) / ...
% 	std(mean(mea.hga(mea.Time < 0, :), 2));
% hga = hga(ts(temp));  % high-gamma at discharge times
% mea.dischargeInds = ts(temp(hga > 2));
% mea.dischargeTimes = mea.Time(mea.dischargeInds);  % peaks with high-gamma at least one sd above the mean
% 
% mea.dischargeTimes(mea.dischargeTimes < mea.epochs(2)) = [];  % exclude recruitment discharges
% mea.dischargeTimes(mea.dischargeTimes > (mea.Time(end) - mea.Padding(2))) = [];  % exclude post termination discharges
% 
% mea.IDIs = diff(mea.waveTimes);
% mea.cov = zeros(floor((numel(mea.IDIs) - 10) / 3), 1);
% for i = 1:numel(mea.cov)
% 	ii = 3*(i-1);
% 	interval = min(30, numel(mea.IDIs) - ii);
% 	temp = mea.IDIs(ii+1:ii+interval);
% 	mea.cov(i) = std(temp) / mean(temp);
% end
% 
% [~, mea.epochs(3)] = max(diff(mea.cov));
% mea.epochs(3) = mea.dischargeTimes(3 * (mea.epochs(3) - 1) + 1);
% mea.dischargeTimes = mea.Time(mea.dischargeInds);
