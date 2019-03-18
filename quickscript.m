% filename = m.Name;
% mea = m.Neuroport;
% pat = 'c5'; seizure = 1;

mea = matfile(sprintf('%s_Seizure%d_Neuroport', pat, seizure), ...
	'writable', true);
mea = filter_mea(mea, [], {'mua'; 'lfp'});
mea = mua_events(mea);
mea = mua_firing_rate(mea);
mea = test_for_recruitment(mea, 'fano');
mea = wave_prop(mea);
disp('Saving figures...')
nn = strrep(mea.Name, 'Seizure', '');
print(2, nn, '-dpng')
print(3, [nn '_FF'], '-dpng')
disp('ALL DONE! GO HOME!!!')
% print(4, [nn, '_coh'], '-dpng')
