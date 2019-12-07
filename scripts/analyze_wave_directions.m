% fname = 'SCM/SCM_Seizure0_Neuroport_10_30.mat'; paramfile = '';
% datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
% addpath(datapath);  % ... add the original data path first
% patpath = genpath(pat);  % ... and then add the local patient path on top 
% addpath(patpath);  % ... so that it is searched first

% Load mea data
[fpath, fname, ~] = fileparts(fname);
addpath(fpath);
mea = load(fname);

% Add field params to mea
if isempty(paramfile)
	mea.params = init_mea_params();
else
	run(paramfile);  % this file should add field ''params'' to mea
end

% Create output file
[~, name, ~] = fileparts(fname);
outname = sprintf('%s_wave_prop', name);
outfile = matfile(outname, 'writable', true);

% Compute discharge times
wt = get_discharge_times(mea);
mea.waveTimes = wt; clear wt

%%

disp('Computing wave directions from events ...')
[events, mea] = wave_prop(mea, 'events');
plot_wave_directions(mea, events);
print(gcf, events.Name, '-dpng');
events.params = mea.params;
outfile.events = events;
% 
disp('Computing wave directions from maxdescent ...')
[maxdescent, mea] = wave_prop(mea, 'maxdescent');
plot_wave_directions(mea, maxdescent);
print(gcf, maxdescent.Name, '-dpng');
maxdescent.params = mea.params;
outfile.maxdescent = maxdescent;

% disp('Computing wave directions from rising deviance ...')
% [rising, mea] = wave_prop(mea, 'rising', 'exclude', false);
% plot_wave_directions(mea, rising);
% print(gcf, rising.Name, '-dpng');
% rising.params = mea.params;
% outfile.rising = rising;
% 
% disp('Computing wave directions from falling deviance ...')
% [falling, mea] = wave_prop(mea, 'falling', 'thresh', -Inf, 'exclude', false);
% plot_wave_directions(mea, falling);
% print(gcf, falling.Name, '-dpng');
% falling.params = mea.params
% outfile.falling = falling;
% 
disp('Computing wave directions from delays ...')
T = .2;
band = [0 50];
mea.params.T = T;
mea.params.delay_band = band;
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, delays);
print(gcf, delays.Name, '-dpng')
delays.params = mea.params;
fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
outfile.(fieldname) = delays;

disp('Computing wave directions from delays (again) ...')
T = 1;
band = [1 50];
mea.params.T = T;
mea.params.delay_band = band;
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, delays);
print(gcf, delays.Name, '-dpng')
delays.params = mea.params;
fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
outfile.(fieldname) = delays;

disp('Computing wave directions from delays (again) ...')
T = 1;
band = [1 13];
mea.params.T = T;
mea.params.delay_band = band;
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, delays);
print(gcf, delays.Name, '-dpng')
delays.params = mea.params;
fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
outfile.(fieldname) = delays;

disp('Computing wave directions from delays (again) ...')
T = 10;
band = [1 50];
mea.params.T = T;
mea.params.delay_band = band;
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, delays);
print(gcf, delays.Name, '-dpng')
delays.params = mea.params;
fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
outfile.(fieldname) = delays;

disp('Computing wave directions from delays (T=10) ...')
T = 10;
band = [1 13];
mea.params.T = T;
mea.params.delay_band = band;
[delays, mea] = wave_prop(mea, 'delays');
plot_wave_directions(mea, delays);
print(gcf, delays.Name, '-dpng')
delays.params = mea.params;
fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
outfile.(fieldname) = delays;


disp('Done.')
disp(outname)

