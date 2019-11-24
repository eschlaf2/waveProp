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

% Get parameters and display them (for remote run tracking)
if ~exist('params', 'var'); params = init_scm_params(); end
disp(params.meta)

create_directory(params);
run_simulation(params);
convert_to_mea(params.meta);

fname = sprintf('%s/SCM/SCM_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, params.meta.sim_num, params.meta.padding);
paramfile = 'tempscripts/scm_mea_options.m';
analyze_wave_directions;

%% Subroutines
function create_directory(params)

	PM = params.meta;
	if ~PM.save, return, end
	
	[base_path, ~, ~] = fileparts(PM.basename);
	if ~exist(base_path, 'dir'), mkdir(base_path), end
	save(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'params');
end

function run_simulation(params)

PM = params.meta;

% Extract variables from meta-parameters
BASENAME = PM.basename;
DURATION = PM.duration;
PADDING = PM.padding;
SAVE = PM.save;
SIM_NUM = PM.sim_num;
T_STEP = PM.t_step;
T0_START = PM.t0_start;

% Create stimulus map
if T0_START > 0
	% If not starting a fresh sim, use previously saved <last>
	load(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, T0_START - 1), 'last')
else
	% ... otherwise, start a fresh sim
	last = params.IC;  %Load the initial conditions to start.
end

K = sum(PADDING) + DURATION;
fig = [];
for t0 = PM.t0_start:T_STEP:K-1  % For each step
	tic
	% Update time offset
	params.t0 = t0 - PADDING(1);
	
	% ... show progress, 
	fprintf('Running %d / %d .. ', t0, floor(K-1));  
		
	% ... get appropriate source drive,
	source_drive = set_source_drive(t0, PM);  
	
	% ... run simulation for duration T_STEP,
	[NP, EC, time, last, fig] = ...  
		seizing_cortical_field(source_drive, min(T_STEP, K - t0 - 1), last, fig, params);
	
	% ... save the results of this run,
	if SAVE
		fprintf('Saving .. ')
		fname = checkname(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, t0*T_STEP));
		save(fname, 'NP','EC','time','last');
	end
	toc
	% ... update progress.
	fprintf('Done.\n')  
end

end

function convert_to_mea(PM)
	files = dir(sprintf('%s_%d_*mat', PM.basename, PM.sim_num));
	addpath(files(1).folder);
	load(files(1).name, 'last');
	cmap = bone;
	im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
	movQ(numel(files) - 1) = im2frame(im, cmap);
	movV(numel(files) - 1) = im2frame(im, cmap);
	[qe, ve, tt] = deal(cell(numel(files) - 1 , 1));
	
	for f = files'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		ind = strsplit(f.name, {'_', '.'});
		ind = str2double(ind{end - 1}) + 1;
		disp(ind)
		im = round(rescale(last.Qe) * (length(cmap) - 1)) + 1;
		movQ(ind) = im2frame(im, cmap);
		movV(ind) = im2frame(round(rescale(last.Ve) * (length(cmap) - 1)) + 1, cmap);
		qe{ind} = NP.Qe;
		ve{ind} = NP.Ve;
		tt{ind} = time;
	end
	
	ve_mat = -cat(1, ve{:});
	qe_mat = cat(1, qe{:});
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, PM.subsample);
	dt = 1 / sample_rate;
	nt = size(ve_mat, 1);
	inds = interp1(time, 1:nt, time(1):dt:time(end), 'nearest');
	time =@() time(1):dt:time(end);
	ve_mat = ve_mat(inds, :, :);
	qe_mat = qe_mat(inds, :, :);
	
	
	mea = create_mea( ...
		ve_mat, ... 
		'SamplingRate', sample_rate, ... 
		'Padding', PM.padding, ...
		'Name', ['SCM Seizure ' num2str(PM.sim_num)], ...
		'Time', time, ... 
		'Path', sprintf('%s/SCM/SCM_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, PM.sim_num, PM.padding) ...	 
		);
	mea.firingRate = reshape(qe_mat, size(mea.Data));
	mea.event_inds = rate2events(mea);
	mea.event_mat_size = size(mea.Data);
	fprintf('Saving %s ... ', mea.Path);
	save(mea.Path, '-struct', 'mea');
	m = matfile(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'Writable', true);
	m.Qe_movie = movQ;
	m.Ve_movie = movV;
	
	fprintf('Done.\n')
	
	fprintf('Done.\n')
	
end

%% Helpers
function event_inds = rate2events(mea)
	lambda = rescale(single(mea.firingRate), 0, 500);  % range is based on MG49_43
	X = rand(size(lambda));
	events = X > exp(-lambda / mea.SamplingRate);
	event_inds = find(events);
end

function [source_drive] = set_source_drive(t, params)

	if t < params.padding(1)  % preseizure
		source_drive = 0;
	elseif t >= params.padding(1) && t < (params.padding(1) + params.duration)  % seizure
		source_drive = params.source_drive; 
	else  % postseizure
		source_drive = params.post_ictal_source_drive;
	end

end

