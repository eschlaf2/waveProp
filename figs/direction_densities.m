% S = 207;
mea = load(sprintf('SCM/SCM_Seizure%d_Neuroport_10_10.mat', S));
fits = load(sprintf('SCM_Seizure%d_Neuroport_10_10_wave_prop.mat', S));

close([1 10])
clear hist_fits;

mea.params = fits.maxdescent.params;
% mea.params = init_mea_params();
time = mea.Time();

fr = mua_firing_rate(mea);
% load('firing_rates/MG49_Seizure43.mat', 'fr')
% mea.Time = downsample(time, 30);

f = {'delays_T10_fband1_13', 'events', 'maxdescent'};
figure(2);
for ii = 1:3; subplot(3,1,ii); plot_dir_simple(mea, fits.(f{ii})); title(rename_metrics(f{ii})); end
% figure(10); plot_dir_simple(mea, fits.events);

compute_times = fits.events.computeTimes;
time_inds = compute_times >= (time(end) - mea.Padding(2) - 30) * 1e3;
% time_inds = compute_times > 0;

 
hist_fits.d10 = fits.delays_T10_fband1_13.Z(time_inds & (fits.delays_T10_fband1_13.p(:) < .05));
% hist_fits.d = fits.delays_T01_fband1_50.Z(time_inds & (fits.delays_T01_fband1_50.p(:) < .05));
hist_fits.e = fits.events.Z(time_inds & (fits.events.p(:) < .05));
hist_fits.m = fits.maxdescent.Z(time_inds & (fits.maxdescent.p(:) < .05));

figure(1) 
for f = fieldnames(hist_fits)'
    [d, xi] = ksdensity(hist_fits.(f{:}));
    plot(xi, d/max(d), 'DisplayName', f{:}); hold on;
end
hold off;
figure(1); legend()
M = load(sprintf('SCM/SCM/SCM_%d_info.mat', S));
figure(1); movie(M.Qe_movie)


