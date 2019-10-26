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
IC = params.IC;
padding = params.padding;
SAVE = params.SAVE;
sim_num = params.sim_num;
t_step = params.t_step;

% Create stimulus map
[map,state] = make_map(params,-padding(1),0);           % Then, make the source map.
last = IC;           %Load the initial conditions to start.

K = sum(padding) + duration;  
fig = [];
for t0 = 0:t_step:K-1
	fprintf('Running %d / %d .. ', t0, K-1);
	[source_drive, map, state] = ...
		set_source_drive(t0, last, params, map, state);

	[NP, EC, time, last, fig] = ...
		seizing_cortical_field(source_drive, map, t_step, last, fig, params);
	time = time - padding(1) + t0;
	
	% Save the results of this run.
	if SAVE
		fprintf('Saving .. ')
		fname = sprintf('%s_%d_%03d.mat', basename, sim_num, t0*t_step);
		save(fname, 'NP','EC','time','last');
	end
	
	fprintf('Done.\n')
end

end

%% Convert to mea
function convert_to_mea(params)
	files = dir(sprintf('%s_%d_*mat', params.basename, params.sim_num));
	addpath(files(1).folder);
	fig = figure();
	ih = imagesc(params.IC.dVe);
	th = title('0');
	clear mov
	mov(numel(files) - 1) = getframe(fig);
	[data, tt] = deal(cell(numel(files) - 1 , 1));
	
	for f = files'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		ind = strsplit(f.name, {'_', '.'});
		ind = str2double(ind{end - 1}) + 1;
		disp(ind)
		set(ih, 'cdata', last.dVe);
		set(th, 'string', num2str(ind)); 
		drawnow
		mov(ind) = getframe(fig);
		data{ind} = NP.Ve;
		tt{ind} = time;
	end
	mov(cellfun(@isempty, {mov.cdata})) = [];
	data(cellfun(@isempty, data)) = [];
	tt(cellfun(@isempty, tt)) = [];
	
	temp = cat(1, data{:});
	mea.Data = reshape(temp)
	
	
	
end

%% Sub routines
function [source_drive, map, state] = set_source_drive(t0, last, params, map, state)

if t0 < params.padding(1)  % preseizure
	source_drive = mean(last.dVe(:));
elseif t0 >= params.padding(1) && t0 < (params.padding(1) + params.duration)  % seizure
	source_drive = params.ictal_source_drive; 
	if strcmp(params.map_type, 'ictal_wavefront')      %... when appropriate, update ictal wavefront location.
		[map,state] = make_map(params,t0,state);
	end
else  % postseizure
	source_drive = params.post_ictal_source_drive;
end
end
