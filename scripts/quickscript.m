% pat = 'c5'; seizure = 3;
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
if ~exist(fname, 'file')
	disp('Creating epoch file ...')
	create_epoch(pat, seizure, 'padding', [10 10]);
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
outfile = matfile([name '_wave_prop_all_waves'], 'writable', true);

% disp('Computing wave directions from events ...')
% [events, mea] = wave_prop(mea, 'events');
% plot_wave_directions(mea, events);
% print(gcf, events.Name, '-dpng');
% outfile.events = events;
% 
% disp('Computing wave directions from maxdescent ...')
% [maxdescent, mea] = wave_prop(mea, 'maxdescent');
% plot_wave_directions(mea, maxdescent);
% print(gcf, maxdescent.Name, '-dpng');
% outfile.maxdescent = maxdescent;

% disp('Computing wave directions from rising deviance ...')
% [rising, mea] = wave_prop(mea, 'rising');
% plot_wave_directions(mea, rising);
% print(gcf, rising.Name, '-dpng');
% outfile.rising = rising;

disp('Computing wave directions from falling deviance ...')
mea = exclude_channels(mea);
[~, mea] = get_discharge_times(mea, 'method', 2);
[falling, mea] = wave_prop(mea, 'falling', 'thresh', -Inf, 'exclude', false);
plot_wave_directions(mea, falling);
print(gcf, falling.Name, '-dpng');
outfile.falling = falling;

% disp('Computing wave directions from delays ...')
% [delays, mea] = wave_prop(mea, 'delays');
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% outfile.delays = delays;

disp('Done.')

