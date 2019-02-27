filename = m.Name;
mea = m.Neuroport;
mea = filter_mea_ecog(mea, [], [], {'mua'});
mea = filter_mea_ecog(mea, [], [], {'lfp'});
mea = mua_events(mea);
mea = mua_firing_rate(mea);
mea = test_for_recruitment(mea, 'fano');
mea = wave_prop(mea);
disp('Saving figures')
nn = strrep(m.Name, 'Seizure', '');
print(2, nn, '-dpng')
print(3, [nn '_FF'], '-dpng')
print(4, [nn, '_coh'], '-dpng')
