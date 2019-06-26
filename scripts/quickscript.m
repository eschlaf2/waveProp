pat = 'c7'; seizure = 1;
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first
computetimesmethod = 1;
showplots = false;

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
outfile = matfile([name '_cohgram'], ...
	'writable', true);
% mea = exclude_channels(mea);
[~, mea] = filter_mea(mea, 'lfp');
% [~, mea] = get_discharge_times(mea, 'method', computetimesmethod);
% mea.Time = mea.Time();

T = 10;  % Window (s)
STEP = 1;  % Step (s)
THRESH = 5e-2;  % significance threshold
TW = 20;  % bandwidth (Hz)
FS = floor(mea.SamplingRate / mea.skipfactor);  % sampling frequency (Hz)
% FPASS = [0 100];  % Frequencies of interest

data1 = mea.lfp;
data2 = mea.lfp;
movingwin = [10 1];  % [window step] seconds
params.err = [1 THRESH];
params.Fs = FS;
% params.fpass
params.tapers = [TW 2*TW-1];  

[C, phi, S12, S1, S2, t, f, confC, phistd] = ...
    cohgramc(data1, data2, movingwin, params);

outfile.C = C;
outfile.phi = phi;
outfile.S12 = S12;
outfile.S1 = S1;
outfile.S2 = S2;
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
outfile.phistd = phistd;


% disp('Computing wave directions from events ...')
% [events, mea] = wave_prop(mea, 'events', ...
% 	'exclude', false, 'showplots', showplots);
% plot_wave_directions(mea, events);
% print(gcf, events.Name, '-dpng');
% outfile.events = events;
% 
% disp('Computing wave directions from maxdescent ...')
% [maxdescent, mea] = wave_prop(mea, 'maxdescent', ...
% 	'exclude', false, 'showplots', showplots);
% plot_wave_directions(mea, maxdescent);
% print(gcf, maxdescent.Name, '-dpng');
% outfile.maxdescent = maxdescent;
% 
% % disp('Computing wave directions from rising deviance ...')
% % [rising, mea] = wave_prop(mea, 'rising', 'exclude', false);
% % plot_wave_directions(mea, rising);
% % print(gcf, rising.Name, '-dpng');
% % outfile.rising = rising;
% % 
% % disp('Computing wave directions from falling deviance ...')
% % [falling, mea] = wave_prop(mea, 'falling', 'thresh', -Inf, 'exclude', false);
% % plot_wave_directions(mea, falling);
% % print(gcf, falling.Name, '-dpng');
% % outfile.falling = falling;
% 
% disp('Computing wave directions from delays ...')
% [delays, mea] = wave_prop(mea, 'delays', ...
% 	'exclude', false, 'showplots', showplots);
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% outfile.delays = delays;

disp('Done.')

