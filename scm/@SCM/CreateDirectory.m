function CreateDirectory(params)

	PM = params.meta;
	if ~PM.save, return, end
	
	[base_path, ~, ~] = fileparts(PM.basename);
	if ~exist(base_path, 'dir'), mkdir(base_path), end
	save(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'params');
end
