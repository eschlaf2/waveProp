% S = 207;
% mea = load(sprintf('SCM/SCM_Seizure%d_Neuroport_10_10.mat', S));
% fits = load(sprintf('SCM_Seizure%d_Neuroport_10_10_wave_prop.mat', S));

% pat = 'MG49'; 
% seizure = 43;

mea = load(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', pat, pat, seizure));
fits = load(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pat, seizure));
mea.params = fits.maxdescent.params;

try close([1 2]); catch ME; end
clear hist_fits;

time = mea.Time();

if ~isfield(mea, 'firingRate'), [~, mea] = mua_firing_rate(mea); end

f = {'delays_T10_fband1_13', 'events', 'maxdescent', 'delays_T01_fband1_13'};

% f = {'delays_T01_fband1_13', 'delays_T10_fband1_50', 'delays_T01_fband1_50'};
figure(2);
for ii = 1:numel(f); subplot(numel(f),1,ii); plot_dir_simple(mea, fits.(f{ii})); title(rename_metrics(f{ii})); end

%%
compute_times = fits.events.computeTimes;
% time_inds = compute_times >= (time(end) - mea.Padding(2) - 30) * 1e3;
time_inds = compute_times > 0;

for ii = 1:numel(f)
    mask = time_inds & isfinite(fits.(f{ii}).p(:)) & (fits.(f{ii}).p(:) < .05);
    ct.(rename_metrics(f{ii})) = compute_times(mask);
    hist_fits.(rename_metrics(f{ii})) = fits.(f{ii}).Z(mask);
end
 
% hist_fits.d10 = fits.delays_T10_fband1_13.Z(time_inds & (fits.delays_T10_fband1_13.p(:) < .05));
% % hist_fits.d = fits.delays_T01_fband1_50.Z(time_inds & (fits.delays_T01_fband1_50.p(:) < .05));
% hist_fits.e = fits.events.Z(time_inds & (fits.events.p(:) < .05));
% hist_fits.m = fits.maxdescent.Z(time_inds & (fits.maxdescent.p(:) < .05));
% 
rotateby = angle(nansum(exp(hist_fits.E * 1j)));
gridx1 = min(compute_times(time_inds))/1e3:.1:max(compute_times(time_inds))/1e3;
gridx2 = linspace(-3*pi, 3*pi, 301);
[x1, x2] = meshgrid(gridx1, gridx2);
x1 = x1(:); x2 = x2(:);

figure(1) 
for f = fieldnames(hist_fits)'
% 	temp = [hist_fits.(f{:}), hist_fits.(f{:}) + 2*pi, hist_fits.(f{:}) - 2*pi]; 
	temp = angle(exp(1j * (hist_fits.(f{:}) - rotateby)));
	temp = [temp(:) - 2*pi; temp(:); temp(:) + 2*pi];
% 	ksdensity([repmat(ct.(f{:}), 3, 1) temp], [x1 x2], 'bandwidth', [5 .15*pi/3]);
	
   [d, xi, bw] = ksdensity(temp, gridx2, 'bandwidth', .15*pi/3);
    plot(xi, d/max(d), 'DisplayName', f{:}); hold on;
	xlim([-pi pi])
end
hold off;
figure(1); legend()
% M = load(sprintf('SCM/SCM/SCM_%d_info.mat', S));
% figure(1); movie(M.Qe_movie)


