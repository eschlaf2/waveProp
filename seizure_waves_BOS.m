function [mea] = seizure_waves_BOS(mea, PLOT)

%% Parse input
if ~exist('PLOT', 'var')
	PLOT = true;
end

% data = mea.lfp;
position = [mea.X mea.Y];

%% Compute coherence and phase
% Note: the chronux toolbox must be on the path

disp('Computing coherence and phase')

if strcmp(which('chronux'), '') % required toolbox
	error('Chronux toolbox not found. Add the toolbox to your matlab path using `addpath(genpath(path_to_chronux))`.');
end

BAND = [1 13];                  % Select a frequency range to analyze
T = 2;							% Length of recording (s)
W = 2;                          % Bandwidth
ntapers = 2*(T * W) - 1;        % Choose the # of tapers.
params.tapers = [T * W, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = mea.SamplingRate / 30;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.

[coh, phi, freq, coh_conf] = compute_coherence(mea.lfp(1:30:(2*3e4), :), params);
% save(['../data/example_' data_set '_coherence'], 'coh', 'phi', 'freq', 'coh_conf');

if PLOT
	figure;
	plot(freq, squeeze(coh(1,8,:)));
	xlabel('Frequency (Hz)');
	ylabel('Coherence');

	figure;
	subplot(1,2,1)
	imagesc(squeeze(coh(:,:, find(freq >= 6, 1))));
	c = colorbar;
	c.Label.String = 'Coherence at 6 Hz';
	xlabel('Electrodes');
	ylabel('Electrodes');
	subplot(1,2,2)
	imagesc(squeeze(phi(:,:, find(freq >= 6, 1))));
	c = colorbar;
	c.Label.String = 'Phase at 6 Hz';
	xlabel('Electrodes');
	ylabel('Electrodes');
end

%% Estimate delays between electrodes

disp('Estimating delays between electrodes')

[delay, delay_ci_lo, delay_ci_up] = compute_delay(coh, coh_conf, phi, freq);

%%
if PLOT
	figure;
	imagesc(1000 * delay);
	c = colorbar;
	c.Label.String = 'Delay (ms)';
	xlabel('Electrodes');
	ylabel('Electrodes');

	% find center electrode
	[~, center] = min((position(:,1) - mean(position(:,1))).^2 + (position(:,2) - mean(position(:,2))).^2);

	figure;
	scatter(position(:,1), position(:,2), 400, 1000 * delay(center,:), 'filled');
	hold on
	plot(position(center,1), position(center,2), 'X');
	c = colorbar;
	c.Label.String = 'Delay to center electrode X (ms)';
	xlabel('Electrode position (mm)');
	ylabel('Electrodes position (mm)');
end


%% Estimate parameters of the waves
[src_dir, speed, ci_dir, ci_sp] = estimate_wave(delay, position, 'plot');
