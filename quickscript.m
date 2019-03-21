% filename = m.Name;
% mea = m.Neuroport;
% pat = 'c5'; seizure = 3;
patpath = genpath(pat);
addpath(patpath);

try
	mea = matfile(sprintf('%s_Seizure%d_Neuroport', pat, seizure), ...
		'writable', true);
catch ME
	m = matfile(sprintf('%s_Seizure%d', pat, seizure), ...
		'writable', true);
	mea = m.Neuroport;
	
end
disp('Computing event times ...')
mua_events(mea);
disp('Computing firing rate ...')
mua_firing_rate(mea);
% mea = test_for_recruitment(mea, 'fano');
disp('Computing wave directions ...')
wave_prop(mea, 'nyc');
wave_prop(mea, 'bos');
% figs_nyc = plot_wave_directions(mea, 'nyc');
disp('Done.')

% disp('Saving figures...')
% nn = strrep(mea.Name, 'Seizure', '');
% print(2, nn, '-dpng')
% print(3, [nn '_FF'], '-dpng')

% disp('ALL DONE! GO HOME!!!')
% print(4, [nn, '_coh'], '-dpng')

rmpath(patpath);
