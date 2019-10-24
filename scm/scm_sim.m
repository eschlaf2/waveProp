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

if exist('options', 'var'), params = init_scm_params(options{:}), else params = init_scm_params(), end
struct2var(params);
[base_path, ~, ~] = fileparts(basename);
if ~exist(base_path, 'dir'), mkdir(base_path), end

%%%% Define the source map. -----------------------------------------------
[map,state] = make_map(params,-padding(1),0);           % Then, make the source map.

%%%% Run the simulation. --------------------------------------------------
last = IC;           %Load the initial conditions to start.

if strcmp(map_type, 'ictal_wavefront')          %When using the ictal wavefront map,     
	last.dVe = zeros(size(last.dVe)) - 1 ;      %... adjust offset to resting potential of excitatory population (mV)
end

K = sum(padding) + duration;  
for t0 = 0:t_step:K
	if t0 <= padding(1)  % preseizure
		source_drive = mean(last.dVe(:));
	elseif t0 > padding(1) && t0 < padding(1) + duration  % seizure
		source_drive = 3; 
		if strcmp(map_type, 'ictal_wavefront')      %... when appropriate, update ictal wavefront location.
			[map,state] = make_map(params,t0,state);
		end
	else  % postseizure
		source_drive = 1.5;
	end

	[NP, EC, time, last] = seizing_cortical_field(source_drive, map, t_step, last, [], params);
	time = time - padding(1) + t0;
	%... save the results of this run.
	save(sprintf('%s_%d_%03d.mat', basename, sim_num, t0*T0), ...
		'NP','EC','time','last')

end




