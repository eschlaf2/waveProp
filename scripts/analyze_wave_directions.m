% pat = 'c7'; seizure = 1;
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first
computetimesmethod = 1;
T = 2;
band = [1 50];
showplots = false;

fname = sprintf('%s_Seizure%d_Neuroport_60_60.mat', pat, seizure);
if ~exist(fname, 'file')
	disp('Creating epoch file ...')
	create_epoch(pat, seizure, 'padding', [60 60]);
end

mea = load(fname);
if ~isfield(mea, 'Path')
	mea.Path = which(fname);
	save(mea.Path, '-v7.3', '-struct', 'mea')
end
if ~isfield(mea, 'BadChannels')
	mea.BadChannels = [];
	save(which(fname), '-v7.3', '-struct', 'mea');
end
[~, name, ~] = fileparts(fname);

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
outname = sprintf('%s_wave_prop_%d', name, computetimesmethod);
outfile = matfile(outname, 'writable', true);
mea = exclude_channels(mea);
[~, mea] = get_discharge_times(mea, 'method', computetimesmethod);

%%

disp('Computing wave directions from events ...')
[events, mea] = wave_prop(mea, 'events', ...
	'exclude', false, 'showplots', showplots);
plot_wave_directions(mea, events);
print(gcf, events.Name, '-dpng');
outfile.events = events;

% disp('Computing wave directions from maxdescent ...')
% [maxdescent, mea] = wave_prop(mea, 'maxdescent', ...
% 	'exclude', false, 'showplots', showplots);
% plot_wave_directions(mea, maxdescent);
% print(gcf, maxdescent.Name, '-dpng');
% outfile.maxdescent = maxdescent;

% disp('Computing wave directions from rising deviance ...')
% [rising, mea] = wave_prop(mea, 'rising', 'exclude', false);
% plot_wave_directions(mea, rising);
% print(gcf, rising.Name, '-dpng');
% outfile.rising = rising;
% 
% disp('Computing wave directions from falling deviance ...')
% [falling, mea] = wave_prop(mea, 'falling', 'thresh', -Inf, 'exclude', false);
% plot_wave_directions(mea, falling);
% print(gcf, falling.Name, '-dpng');
% outfile.falling = falling;

% disp('Computing wave directions from delays ...')
% [delays, mea] = wave_prop(mea, 'delays', ...
% 	'exclude', false, 'showplots', showplots, 'T', T, 'band', band);
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% outfile.(sprintf('delays_T%02d_fband%d_%d', T, band)) = delays;

disp('Done.')

