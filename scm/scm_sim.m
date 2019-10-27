%Run the Seizing Cortical Model with a chosen source distribution map.
%
%CALLS:
%
% . seizing_cortical_field.m
% . make_map.m
% . init_scm_params.m
%
%OUTPUTS:
%
% . The output of function "seizing_cortical_field.m" is saved after
%   each 1 s simulation.  The output is saved in a .mat file 
%   with name "seizing_cortical_field_" MAP_TYPE "k_X.mat" where "X" is an
%   integer specifying the time step, and MAP_TYPE is either "fixed_point_source"
%   or "ictal_wavefront".  The file is saved under basename,
%   specified in init_scm_params.m.
%  
% . visualize_scm_excitatory_population:
%   Plots an image of the excitatory population
%   activity every 2 s. To do so, the code reads in the .mat files
%   generated in this cell and saved according to basename.
%
% . save_scm_ecog_data: Select a 10 s interval of simulated ECOG data, and convert it
%   to a format appropriate for wave analysis. **UPDATE THIS TO CONVERT TO
%   MEA**
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seizing Cortical Model (a modification of the Waikato Cortical Model)
% Copyright (c) 2016 M. A. Kramer
% Department of Mathematics and Statistics, Boston University, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('options', 'var'); options = {}; end

params = init_scm_params(options{:});
disp(params)

[base_path, ~, ~] = fileparts(params.basename);
if ~exist(base_path, 'dir'), mkdir(base_path), end
if params.SAVE, save(sprintf('%s_%d_info', params.basename, params.sim_num), 'params'); end

run_simulation(params)
convert_to_mea(params)

%% Run Simulation
function run_simulation(params)

% Extract parameters
basename = params.basename;
duration = params.duration;
padding = params.padding;
SAVE = params.SAVE;
sim_num = params.sim_num;
t_step = params.t_step;
t0_start = params.t0_start;

% Create stimulus map
if t0_start > 0
	load(sprintf('%s_%d_%03d', basename, sim_num, t0_start - 1), 'last')
else
	last = params.IC;  %Load the initial conditions to start.
end

K = sum(padding) + duration;  
fig = [];
for t0 = t0_start:t_step:K-1
	fprintf('Running %d / %d .. ', t0, K-1);
	last.t0 = t0;
	source_drive = set_source_drive(t0, last, params);

	[NP, EC, time, last, fig] = ...
		seizing_cortical_field(source_drive, t_step, last, fig, params);
	time = time - padding(1);
	
	% Save the results of this run.
	if SAVE
		fprintf('Saving .. ')
		fname = checkname(sprintf('%s_%d_%03d', basename, sim_num, t0*t_step));
		save(fname, 'NP','EC','time','last');
	end
	
	fprintf('Done.\n')
end

end

%% Convert to mea

function convert_to_mea(params)
	files = dir(sprintf('%s_%d_*mat', params.basename, params.sim_num));
	addpath(files(1).folder);
	cmap = bone;
	im = round(rescale(params.IC.Ve) * (length(cmap) - 1)) + 1;
	mov(numel(files) - 1) = im2frame(im, cmap);
	[data, tt] = deal(cell(numel(files) - 1 , 1));
	
	for f = files'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		ind = strsplit(f.name, {'_', '.'});
		ind = str2double(ind{end - 1}) + 1;
		disp(ind)
		im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
		mov(ind) = im2frame(im, cmap);
		data{ind} = NP.Qe;
		tt{ind} = time;
	end
	mov(cellfun(@isempty, {mov.cdata})) = [];
	data(cellfun(@isempty, data)) = [];
	tt(cellfun(@isempty, tt)) = [];
	
	data_mat = cat(1, data{:});
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, params.subsample_rate);
	dt = 1 / sample_rate;
	nt = size(data_mat, 1);
	inds = interp1(time, 1:nt, time(1):dt:time(end), 'nearest');
	time =@() time(1):dt:time(end);
	data_mat = data_mat(inds, :, :);
	
	
	mea = create_mea( ...
		data_mat, ... 
		'SamplingRate', sample_rate, ... 
		'Padding', params.padding, ...
		'Name', ['SCM Seizure ' num2str(params.sim_num)], ...
		'Time', time, ... 
		'Path', sprintf('%s/SCM/SCM_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, params.sim_num, params.padding) ...	 
		);
	
	save(mea.Path, '-struct', 'mea');
	m = matfile(sprintf('%s_%d_info', params.basename, params.sim_num), 'Writable', true);
	m.Ve_movie = mov;
	
end

%% Sub routines
function [source_drive] = set_source_drive(t, last, params)

if t < params.padding(1)  % preseizure
	source_drive = mean(last.dVe(:));
elseif t >= params.padding(1) && t < (params.padding(1) + params.duration)  % seizure
	source_drive = params.ictal_source_drive; 
else  % postseizure
	source_drive = params.post_ictal_source_drive;
end

end

