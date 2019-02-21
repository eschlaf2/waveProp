filename = m.Name;
mea = m.Neuroport;
mea = filter_mea_ecog(mea, [], [], {'mua'});
mea = filter_mea_ecog(mea, [], [], {'lfp'});
mea = mua_events(mea);
mea = mua_firing_rate(mea);