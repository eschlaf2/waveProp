% default_params
function res = init_scm_params(varargin)

%% parsing functions
validate = @(x, all) any(validatestring(x, all));  % validate strings

%% Meta
% Save and visualize

G = inputParser;
p = @(varargin) addParameter(G, varargin{:});  % convenience function
G.KeepUnmatched = true;

p('basename', 'SCM/SCM/SCM');
p('sim_num', 0);
p('SAVE', true);  % Save output
p('visualize_results', false);  %Set this variable to true to create plots during simulation.
p('visualization_rate', 10);  % Show this many frames per second
p('t_step', 1);  % Simulate t_step second intervals
p('t0_start', 0);  % Continue from previous sim (IC will use last from file number t0_start-1)

parse(G, varargin{:});
res.meta = G.Results;

%% Seizure definition
G = inputParser;
p = @(varargin) addParameter(G, varargin{:});  % convenience function
G.KeepUnmatched = true;

p('duration', 90);  % Seizure duration
p('padding', [10 30], @(x) isnumeric(x) & numel(x) == 2);  % Padding before and after seizure

% Model parameters
allMapTypes = {'ictal_wavefront', 'fixed_point_source'};
p('map_type', 'ictal_wavefront', @(x) validate(x, allMapTypes));  % Ictal source 
p('stim_center', [0 0], @(x) numel(x) == 2);
p('grid_size', [100 100], @(x) ~mod(x, 2) && numel(x) == 2);  % size of grid to simulate (must be even)
p('dt', 4e-4);
p('Laplacian', [0 1 0; 1 -4 1; 0 1 0]);
p('ictal_source_drive', 3);
p('post_ictal_source_drive', 1.5);
p('subsample_rate', Inf);  % Allow subsampling when making mea
p('spatial_resolution', .3);  % (mm) 
p('expansion_rate', 3);  % every 3 seconds expand the wavefront
p('IC', {});  % Initial conditions of the sim (enter options as name-value pairs)
p('HL', {});  % Constants that define the locations of steady states

parse(G, varargin{:});
model = G.Results;
model = clean_model(model);


%% Noise
G = inputParser; G.KeepUnmatched = true;

p(G, 'noise', 0.5);  % Noise level
p(G, 'noise_sf', []);    % noise scale-factor
p(G, 'noise_sc', 0.2);             % subcortical noise

G.parse(varargin{:});
noise = G.Results;
if isempty(noise.noise_sf); noise.noise_sf = 0.2*20*noise; end

%% Initial conditions
% These are the same as in Waikato-Kramer except that dVe defaults to -1
% instead of 1 (original IC in <map_type='ictal_wavefront'> scenario).
% Additionally, resizing is done via bootstrapping in the case of a grid
% size different from 100x100.

G = inputParser; G.CaseSensitive = true;
p = @(varargin) addParameter(G, varargin{:});  % convenience function

IC = load('default_scm_ICs.mat');
IC = resize(IC, res.scm);

p('D11', IC.D11)
p('D22', IC.D22)
p('F_ee', IC.F_ee)
p('F_ei', IC.F_ei)
p('F_ie', IC.F_ie)
p('F_ii', IC.F_ii)
p('K', IC.K)
p('Phi_ee', IC.Phi_ee)
p('Phi_ei', IC.Phi_ei)
p('Phi_ie', IC.Phi_ie)
p('Phi_ii', IC.Phi_ii)
p('Qe', IC.Qe)
p('Qi', IC.Qi)
p('Ve', IC.Ve)
p('Vi', IC.Vi)
p('dVe', IC.dVe)  % Excitability of excitatory population
p('dVi', IC.dVi)
p('phi2_ee', IC.phi2_ee)
p('phi2_ei', IC.phi2_ei)
p('phi_ee', IC.phi_ee)
p('phi_ei', IC.phi_ei)

[map, state] = update_map(seizure);
p('map', map);
p('state', state);

parse(G, model.IC{:});

IC = G.Results;

%% Time constants
G = inputParser; G.KeepUnmatched = true;
p = @(G, varargin) addParameter(G, varargin{:});

p(G, 'tau_dD', 200);  %inhibitory gap junction time-constant (/s).
p(G, 'tau_dVe', 250);  %excitatory population resting voltage time-constant (/s).
p(G, 'tau_dVi', 250);  %inhibitory population resting voltage time-constant (/s).

parse(G, varargin{:})
time_constants = G.Results;

%% Steady state constants
G = inputParser; G.CaseSensitive = true;

HL = SCM_init_globs;  % See this function for details
for f = fieldnames(HL)', p(G, f{:}, HL.(f{:})); end
G.parse(model.HL{:});

HL = G.Results;

%% Potassium parameters
G = inputParser; G.KeepUnmatched = true;

p(G, 'tau_K', 200);    %time-constant (/s).
p(G, 'k_decay', 0.1);  %decay rate (/s).
p(G, 'kD', 1);         %diffusion coefficient (cm^2/s).
p(G, 'KtoVe', 10);     %impact on excitatory population resting voltage.
p(G, 'KtoVi', 10);     %impact on inhibitory population resting voltage.
p(G, 'KtoD', -50);    %impact on inhibitory gap junction strength.
p(G, 'kR', 0.15);   %scale reaction term. 

parse(G, varargin{:})
potassium = G.Results;

%% Electrode parameters
G = inputParser; G.KeepUnmatched = true;

% Microscale
p(G, 'xNP', []);
p(G, 'yNP', []);
p(G, 'dimsNP', [10 10]);

% Macroscale
p(G, 'xEC', []);
p(G, 'yEC', []); 
p(G, 'scaleEC', 4);
p(G, 'dimsEC', [3 3]);

parse(G, varargin{:});
electrodes = clean_electrodes(G.Results, model);


end

function s = clean_electrodes(s, model)

if isempty(s.xNP), s.xNP = model.Nx / 2; end
if isempty(s.yNP), s.yNP = model.Ny / 2; end
if isempty(s.xEC), s.xEC = model.Nx / 2; end
if isempty(s.yEC), s.yEC = model.Ny / 2; end

end

function s = clean_model(s, HL, IC, model)

% set stimulus center if unset
stim_center = s.stim_center;
if any(stim_center<=0), stim_center = round(s.grid_size * .4); end
s.stim_center = round(stim_center);

% set time resolution
D2 = IC.D22 * model.spatial_resolution.^2;
if s.dt > 2e-4 && (any(D2 >= 0.87) || HL.v > 140)
	s.dt = 2e-4;
	warning('High D2 or HL.v; setting dt to 2e-4.');
end


end

function IC = resize(IC, res)
% Resize IC fields by bootstrapping from original values. The same subset
% of indices is used for each field. i.e. the ICs on the whole are a
% sub(super)-set of the default ICs.

if ~all(res.grid_size == size(IC.D11))
	rng(0)  % for reproducibility
	inds = randi(numel(IC.D11), res.grid_size(1), res.grid_size(2));
	for f = fieldnames(IC)'
		IC.(f{:}) = IC.(f{:})(inds);
	end
end

end