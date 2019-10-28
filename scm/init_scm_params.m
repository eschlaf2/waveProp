% default_params
function res = init_scm_params(varargin)
% All parameters except IC and HL can be adjusted by passing name-value
% pairs. To modify HL and IC, the value passed to either name should be a
% cell containing name-value pairs of specific subparameters.
% 
% Example:
% 	params = init_scm_params( ...
% 		'sim_num', 1, ...
% 		'IC', {'dVe', 1} ...
% 		);

%%%% Adjust values as name-value pairs %%%%
res.meta = parse_meta(varargin); 
	% basename             SAVE                 sim_num
	% t0_start             t_step               visualization_rate
	% visualize            padding              duration
	% subsample            ictal_source_drive   post_ictal_source_drive         

res.model = parse_model(varargin);
	% dt              electrodes      expansion_rate  grid_size       
	% Laplacian       noise           stim_center     spatial_resolution
	% time_constants  SS              IC              K               

	
%%%% Passed as unmatched parameters from parsing of model %%%%
res.electrodes = parse_electrodes(res.model);
	% centerEC        centerNP        dimsEC          dimsNP          scaleEC         

res.time_constants = parse_time(res.model);
	% tau_dD          tau_dVe         tau_dVi         
	
res.K = parse_potassium(res.model);
	% k_decay         kD              kR              KtoD            KtoVe
	% KtoVi           tau_K

res.noise = parse_noise(res.model);
	% noise           noise_sc        noise_sf        

	
%%%% IC (initial conditions) and SS (steady states) %%%%
res.IC = parse_IC(res.model);  
	% D11             D22             F_ee            F_ei            F_ie
	% F_ii            K               Phi_ee          Phi_ei
	% Phi_ie          Phi_ii          Qe              Qi              Ve
	% Vi              dVe             dVi             map
	% phi2_ee         phi2_ei         phi_ee          phi_ei          state


res.SS = parse_SS(res.model);  
	% Lambda          Nee_a           Nee_ab          Nee_b
	% Nee_sc          Nei_a           Nei_ab          Nei_b
	% Nei_sc          Nie_b           Nii_b           Qe_max
	% Qi_max          Ve_rest         Ve_rev          Vi_rest
	% Vi_rev          d_psi_ee        d_psi_ei        d_psi_ie
	% d_psi_ii        gamma_e         gamma_i         phi_ee_sc
	% phi_ei_sc       rho_e           rho_i           sigma_e
	% sigma_i         tau_e           tau_i           theta_e
	% theta_i         v


res.model = check_time_resolution(res.model, res.IC, res.SS);

end

function meta = parse_meta(options)
% Save and visualize
[G, p] = get_parser();

p('basename', 'SCM/SCM/SCM');
p('sim_num', 0);
p('SAVE', true);  % Save output
p('visualize', false);  %Set this variable to true to create plots during simulation.
p('visualization_rate', 10);  % Show this many frames per second
p('t_step', 1);  % Simulate t_step second intervals
p('t0_start', 0);  % Continue from previous sim (IC will use last from file number t0_start-1)
p('duration', 90);  % Seizure duration
p('padding', [10 30], @(x) isnumeric(x) & numel(x) == 2);  % Padding before and after seizure
p('ictal_source_drive', 3);
p('post_ictal_source_drive', 1.5);
p('subsample', Inf);  % Allow subsampling when making mea

parse(G, options{:});
meta = G.Results;
end

function model = parse_model(options)

[G, p] = get_parser();

% Model parameters
allMapTypes = {'ictal_wavefront', 'fixed_point_source'};
p('stim_center', [0 0], @(x) numel(x) == 2);
p('grid_size', [100 100], @(x) ~mod(x, 2) && numel(x) == 2);  % size of grid to simulate (must be even)
p('dt', 2e-4);
p('Laplacian', [0 1 0; 1 -4 1; 0 1 0]);
p('spatial_resolution', .3);  % (mm) 
p('expansion_rate', 3, @(x) x > 0);  % every 3 seconds expand the wavefront; set to Inf for fixed source
p('IC', {});  % Initial conditions of the sim (enter options as name-value pairs)
p('SS', {});  % Constants that define the locations of steady states
p('time_constants', {});  % tau parameters controlling rates
p('K', {});  % Potassium (K-related) parameters
p('noise', {});  % Noise parameters
p('electrodes', {});  % electrode positions


parse(G, options{:});
model = G.Results;
if isempty(model.time_constants), model.time_constants = G.Unmatched; end
if isempty(model.K), model.K = G.Unmatched; end
if isempty(model.electrodes), model.electrodes = G.Unmatched; end
if isempty(model.noise), model.noise = G.Unmatched; end

% set stimulus center if unset
stim_center = model.stim_center;
if any(stim_center<=0), stim_center = round(model.grid_size * .4); end
model.stim_center = round(stim_center);

end

function noise = parse_noise(model)
%% Noise
[G, p] = get_parser();
options = model.noise;

p('noise', 0.5);  % Noise level
p('noise_sf', []);  % noise scale-factor
p('noise_sc', 0.2);  % subcortical noise

if isstruct(options), G.parse(options), else, G.parse(options{:}); end
noise = G.Results;
if isempty(noise.noise_sf); noise.noise_sf = 0.2*20*noise.noise; end

end

function IC = parse_IC(model)
%% Initial conditions
% These are the same as in Waikato-Kramer except that dVe defaults to -1
% instead of 1 (original IC in <map_type='ictal_wavefront'> scenario).
% Additionally, resizing is done via bootstrapping in the case of a grid
% size different from 100x100.

G = inputParser; G.CaseSensitive = true;
p = @(varargin) addParameter(G, varargin{:});  % convenience function

IC = load('default_scm_ICs.mat');

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

[map, state] = update_map(model);
p('map', map);
p('state', state);

parse(G, model.IC{:});

IC = G.Results;
IC = resize(IC, model);

end

function model = check_time_resolution(model, IC, HL)
% set time resolution
D2 = IC.D22 * model.spatial_resolution.^2;
if model.dt > 2e-4 && (any(D2(:) >= 0.87) || HL.v(:) > 140)
	model.dt = 2e-4;
	warning('High D2 or HL.v; setting dt to 2e-4.');
end
end

function time_constants = parse_time(model)

%% Time constants
[G, p] = get_parser();
options = model.time_constants;

p('tau_dD', 200);  %inhibitory gap junction time-constant (/s).
p('tau_dVe', 250);  %excitatory population resting voltage time-constant (/s).
p('tau_dVi', 250);  %inhibitory population resting voltage time-constant (/s).

if isstruct(options), parse(G, options), else, parse(G, options{:}), end
time_constants = G.Results;

end

function SS = parse_SS(model)
% Steady state constants

[G, p] = get_parser(); 
G.CaseSensitive = true;
options = model.SS;

SS = SCM_init_globs;  % See this function for details
for f = fieldnames(SS)', p(f{:}, SS.(f{:})); end

if isstruct(options), G.parse(options), else, G.parse(options{:}); end

SS = G.Results;

end

function potassium = parse_potassium(model)
% Potassium parameters

[G, p] = get_parser();
options = model.K;

p('tau_K', 200);    %time-constant (/s).
p('k_decay', 0.1);  %decay rate (/s).
p('kD', 1);         %diffusion coefficient (cm^2/s).
p('KtoVe', 10);     %impact on excitatory population resting voltage.
p('KtoVi', 10);     %impact on inhibitory population resting voltage.
p('KtoD', -50);    %impact on inhibitory gap junction strength.
p('kR', 0.15);   %scale reaction term. 

if isstruct(options), parse(G, options), else, parse(G, options{:}), end
potassium = G.Results;

end

function electrodes = parse_electrodes(model)

%% Electrode parameters
[G, p] = get_parser();

% Microscale
p('centerNP', []);
p('dimsNP', [10 10]);

% Macroscale
p('centerEC', []);
p('scaleEC', 4);
p('dimsEC', [3 3]);

if isstruct(model.electrodes), parse(G, model.electrodes), else, parse(G, model.electrodes{:}); end
electrodes = clean_electrodes(G.Results, model);

end

function s = clean_electrodes(s, model)
% Make sure electrodes are in grid bounds

NP = s.centerNP;
EC = s.centerEC;
gs = model.grid_size;

if isempty(NP), NP = gs ./ 2; end
if isempty(EC), EC = gs ./ 2; end
% Keep in bounds
NP = min(max(NP, ceil(s.dimsNP / 2)), gs - floor(s.dimsNP / 2));
EC = min(max(EC, ceil(s.dimsEC / 2)), gs - floor(s.dimsEC / 2));

s.centerNP = NP; s.centerEC = EC;
end

function IC = resize(IC, model)
% Resize IC fields by bootstrapping from original values. The same subset
% of indices is used for each field. i.e. the ICs on the whole are a
% sub(super)-set of the default ICs.

if ~all(model.grid_size == size(IC.D11))
	rng(0)  % for reproducibility
	inds = randi(numel(IC.D11), model.grid_size(1), model.grid_size(2));
	for f = fieldnames(IC)'
		IC.(f{:}) = IC.(f{:})(inds);
	end
end

end

function [G, p] = get_parser()
G = inputParser; G.KeepUnmatched = true;
p =@(varargin) addParameter(G, varargin{:});

end