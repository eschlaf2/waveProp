mea = load('SCM/SCM_Seizure16_Neuroport_10_10.mat');
fits = load('SCM_Seizure16_Neuroport_10_10_wave_prop.mat');
% mea = load('FHN/FHN_Seizure0_Neuroport_10_10.mat');
% fits = load('FHN_Seizure0_Neuroport_10_10_wave_prop.mat');
% mea = load('SIM/SIM_Seizure11_Neuroport_10_10.mat');
% fits = load('SIM_Seizure11_Neuroport_10_10_wave_prop.mat');
% mea = load('MG49/MG49_Seizure43_Neuroport_10_10.mat');
% fits = load('MG49_Seizure43_Neuroport_10_10_wave_prop.mat');

% mea.params = fits.maxdescent.params;
mea.params = init_mea_params();

fr = mua_firing_rate(mea);
% load('firing_rates/MG49_Seizure43.mat', 'fr')

time = mea.Time();
compute_times = fits.events.computeTimes;
time_inds = compute_times >= (time(end) - mea.Padding(2) - 30) * 1e3;
% time_inds = compute_times > 0;
 
hist_fits.d10 = fits.delays_T10_fband1_13.Z(time_inds & (fits.delays_T10_fband1_13.p(:) < .05));
hist_fits.d = fits.delays_T01_fband1_50.Z(time_inds & (fits.delays_T01_fband1_50.p(:) < .05));
hist_fits.e = fits.events.Z(time_inds & (fits.events.p(:) < .05));
hist_fits.m = fits.maxdescent.Z(time_inds & (fits.maxdescent.p(:) < .05));
figure(1) 
for f = fieldnames(hist_fits)'
    [d, xi] = ksdensity(hist_fits.(f{:}));
    plot(xi, d/max(d), 'DisplayName', f{:}); hold on;
end
hold off;
legend

mea.firingRate = mean(fr, 2);
figure(10); plot_dir_simple(mea, fits.delays_T10_fband1_13);

