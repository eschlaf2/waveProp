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
if ~exist('params', 'var'); params = SCMParams; end
disp(params.meta)
try rng(params.meta.seed); catch ME; end

create_directory_(params);
run_simulation_(params);
convert_to_mea_(params);

% analyze_wave_directions_(params);

%% Subroutines
function create_directory_(params)

	PM = params.meta;
	if ~PM.save, return, end
	
	[base_path, ~, ~] = fileparts(PM.basename);
	if ~exist(base_path, 'dir'), mkdir(base_path), end
	save(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'params');
end

function run_simulation_(params)

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
if T0_START > -PADDING(1)
	% If not starting a fresh sim, use previously saved <last>
    try
        load(sprintf('%s_%d_%03d', BASENAME, SIM_NUM, T0_START - 1), 'last')
    catch ME
        
    end
else
	% ... otherwise, start a fresh sim
	last = params.IC;  %Load the initial conditions to start.
end

K = PADDING(2) + DURATION;
fig = [];
for t0 = PM.t0_start:T_STEP:K-T_STEP  % For each step
	tic
	% Update time offset
	params.t0 = t0;
	
	% ... show progress, 
	fprintf('Running %d / %d .. ', t0+1, floor(K));  
		
	% ... run simulation for duration T_STEP,
	[NP, EC, time, last, fig] = ...  
		seizing_cortical_field('legacy variable', ...
            min(T_STEP, K - t0), last, fig, params);
	
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

function convert_to_mea_(PM)
    if ~PM.save, return, end
	if ~isfield(PM, 'label'), PM.label = 'SCM'; end
	files = dir(sprintf('%s_%d_*mat', PM.basename, PM.sim_num));
	addpath(files(1).folder);
	load(files(1).name, 'last');
	cmap = bone;
	im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
	movQ(numel(files) - 1) = im2frame(im, cmap);
	movV(numel(files) - 1) = im2frame(im, cmap);
	[qe, ve, tt, file_inds] = deal(cell(numel(files) - 1 , 1));

    file_inds = cellfun(@(f) strsplit(f, {'_', '.'}), {files.name}, 'uni', 0);
    file_inds = cellfun(@(f) str2double(strrep(f{end - 1}, 'M', '-')), file_inds);
    [~, file_order] = sort(file_inds, 'ascend');
    
	ii = 1;
    for f = files(file_order)'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		disp(ii)
		im = round(rescale(last.Qe) * (length(cmap) - 1)) + 1;
		movQ(ii) = im2frame(im, cmap);
		movV(ii) = im2frame(round(rescale(last.Ve) * (length(cmap) - 1)) + 1, cmap);
		qe{ii} = NP.Qe;
		ve{ii} = NP.Ve;
		tt{ii} = time;
        ii = ii + 1;
	end
	
% 	ve_mat = -cat(1, ve{:});
	qe_mat = cat(1, qe{:});
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, PM.subsample);
	dt = 1 / sample_rate;
	nt = size(qe_mat, 1);
	inds = interp1(time, 1:nt, time(1):dt:time(end), 'nearest');
	time =@() time(1):dt:time(end);
% 	ve_mat = ve_mat(inds, :, :);
	qe_mat = qe_mat(inds, :, :);
	
	
	mea = create_mea( ...
		-qe_mat, ... 
		'SamplingRate', sample_rate, ... 
		'Padding', PM.padding, ...
		'Name', [PM.label ' Seizure ' num2str(PM.sim_num)], ...
		'Time', time, ... 
		'Path', sprintf('%s/%s/%s_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, PM.label, PM.label, PM.sim_num, PM.padding) ...	 
		);
	
%   mea = add_noise_(mea, 2);  % Add 3D brownian noise with snr=2; the spectra after this transformation looked similar to recorded seizures - could also use (much) higher snr pink noise
	qe_mat = rescale(single(qe_mat), 0, 25);  % range is based on experimentation
	mea.firingRate = reshape(qe_mat, size(mea.Data));
	mea.event_inds = rate2events_(mea);
	mea.event_mat_size = size(mea.Data);
	mea.params = init_mea_params();
	fprintf('Saving %s ... ', mea.Path);
	save(mea.Path, '-struct', 'mea');
	m = matfile(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'Writable', true);
	m.Qe_movie = movQ;
	m.Ve_movie = movV;
	
	fprintf('Done.\n')
	
	fprintf('Done.\n')
	
end

function analyze_wave_directions_(params)
fname = sprintf('%s/%s/%s_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, params.label, params.label, params.sim_num, params.padding);
paramfile = '';
analyze_wave_directions;
end

%% Helpers

function mea = add_noise_(mea, snr_dB)
	signal = mea.Data;
	brown_noise = randnd(-2, [length(signal), max(mea.Position)]);  % brownian noise: -2; pink noise: -1
	noise = brown_noise(:, mea.locs);
	signal_power = rms(signal);
	noise_power = rms(noise);
	scale_factor = signal_power ./ (10^(log10(snr_dB)/2) .* noise_power);
	noisy_sig = signal + scale_factor .* noise;
	mea.Data = noisy_sig;

end


