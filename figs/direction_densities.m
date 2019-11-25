% S = 207;
% mea = load(sprintf('SCM/SCM_Seizure%d_Neuroport_10_10.mat', S));
% fits = load(sprintf('SCM_Seizure%d_Neuroport_10_10_wave_prop.mat', S));

pat = 'MG49'; 
seizure = 43;

mea = load(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', pat, pat, seizure));
fits = load(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pat, seizure));
mea.params = fits.maxdescent.params;

try close([1 2]); catch ME; end
clear hist_fits;

time = mea.Time();

[~, mea] = mua_firing_rate(mea);

f = {'delays_T10_fband1_13', 'events', 'maxdescent'};
figure(2);
for ii = 1:3; subplot(3,1,ii); plot_dir_simple(mea, fits.(f{ii})); title(rename_metrics(f{ii})); end

%%
compute_times = fits.events.computeTimes;
time_inds = compute_times >= (time(end) - mea.Padding(2) - 30) * 1e3;
% time_inds = compute_times > 0;

 
hist_fits.d10 = fits.delays_T10_fband1_13.Z(time_inds & (fits.delays_T10_fband1_13.p(:) < .05));
% hist_fits.d = fits.delays_T01_fband1_50.Z(time_inds & (fits.delays_T01_fband1_50.p(:) < .05));
hist_fits.e = fits.events.Z(time_inds & (fits.events.p(:) < .05));
hist_fits.m = fits.maxdescent.Z(time_inds & (fits.maxdescent.p(:) < .05));

figure(1) 
for f = fieldnames(hist_fits)'
	temp = [hist_fits.(f{:}), hist_fits.(f{:}) + 2*pi, hist_fits.(f{:}) - 2*pi]; 
% 	temp = hist_fits.(f{:});
    [d, xi] = ksdensity(temp, 'bandwidth', 1.06);  % This bw came from running the function with [-pi pi] support
    plot(xi, d/max(d), 'DisplayName', f{:}); hold on;
	xlim([-pi pi])
end
hold off;
figure(1); legend()
% M = load(sprintf('SCM/SCM/SCM_%d_info.mat', S));
% figure(1); movie(M.Qe_movie)


