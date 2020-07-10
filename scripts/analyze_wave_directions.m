
% Load mea and get times to fit waves
mea = MEA(fname);
outname = [mea.Name '_fits.mat'];
tAll = 0:.01:mea.Time(end) - mea.Padding(2);
tWave = mea.get_wave_times('events');
outname = [mea.Name '_fits.mat'];
fits = struct();
if ~exist(outname, 'file')
	save(outname, '-struct', 'fits');
end
% FR = load(mea.Path, 'Data');
% FR = -single(FR.Data);
% FR = (FR - min(FR)) ./ std(FR);
% [~, event_inds, ~, ~] = findpeaks(FR(:), 'MinPeakHeight', 2); 
% save(mea.Path, 'event_inds', '-append');
% mea.event_inds = event_inds;

% rng(str2double(mea.seizure));
% mea.SamplingRate = mea.SRO;
% Data = mea.add_noise(2);
% mea.Raw = Data; mea.SamplingRate = 1e3;
% save(mea.Path, 'Data', '-append');
% mea = MEA(mea.Path);

S = warning;
warning('off');
tic

fits.M = mea.max_descent(tAll);
fits.Mdt = mea.max_descent(tWave);
fits.E = mea.dir_events(tAll);
fits.Edt = mea.dir_events(tWave);
save(outname, '-struct', 'fits', '-append');
toc
% tic
% fits.D1xwh = mea.delays(tAll, 'halfwin', .5, 'maskby', 'highest', 'fband', [1 50]);
% % fits.D10 = mea.delays(tAll, 'halfwin', 5);
% save(outname, '-struct', 'fits', '-append');
% toc


	
disp('Success!!')
warning(S);


% May need to do D1 and D10 separately as they take a long time
% #!/bin/bash
% for ss in d1 d10; do qsub -N SW_$ss -t 101-200 run_fits_sims.sh ${ss}.m FS; done


%%

if 0  % old
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
% tt = mea.Time();
% wt = 0:.01:tt(end) - 10;
% mea.params.min_peak_prominence = 0;
% mea.params.min_peak_distance = 10;
tWave = get_discharge_times(mea);
mea.waveTimes = tWave; clear wt

%%

disp('Computing wave directions from events ...')
[events, mea] = wave_prop(mea, 'events');
plot_wave_directions(mea, events);
print(gcf, events.Name, '-dpng');
events.params = mea.params;
outfile.E = events;
% 
disp('Computing wave directions from maxdescent ...')
[maxdescent, mea] = wave_prop(mea, 'maxdescent');
plot_wave_directions(mea, maxdescent);
print(gcf, maxdescent.Name, '-dpng');
maxdescent.params = mea.params;
outfile.M = maxdescent;

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
% outfile.D1xw = delays;

disp('Computing wave directions from delays ...')
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
% outfile.D1 = delays;

% disp('Computing wave directions from delays (again) ...')
% T = 1;
% band = [1 50];
% mea.params.T = T;
% mea.params.delay_band = band;
% [delays, mea] = wave_prop(mea, 'delays');
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% delays.params = mea.params;
% fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
% outfile.(fieldname) = delays;
% 
% disp('Computing wave directions from delays (again) ...')
% T = 1;
% band = [1 13];
% mea.params.T = T;
% mea.params.delay_band = band;
% [delays, mea] = wave_prop(mea, 'delays');
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% delays.params = mea.params;
% fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
% outfile.(fieldname) = delays;
% 
% disp('Computing wave directions from delays (again) ...')
% T = 10;
% band = [1 50];
% mea.params.T = T;
% mea.params.delay_band = band;
% [delays, mea] = wave_prop(mea, 'delays');
% plot_wave_directions(mea, delays);
% print(gcf, delays.Name, '-dpng')
% delays.params = mea.params;
% fieldname = checkname(sprintf('delays_T%02g_fband%d_%d', T, band));
% outfile.(fieldname) = delays;

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
outfile.D10 = delays;


disp('Done.')
disp(outname)

end

