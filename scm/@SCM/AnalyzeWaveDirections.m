function AnalyzeWaveDirections(params)
fname = sprintf('%s/%s/%s_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, params.label, params.label, params.sim_num, params.padding);
paramfile = '';
analyze_wave_directions;
end
