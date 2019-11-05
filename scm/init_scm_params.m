function res = init_scm_params(varargin)
% All parameters except IC and SS can be adjusted by passing name-value
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
		% basename             save                 sim_num
		% t0_start             t_step               visualization_rate
		% padding              duration
		% subsample            source_drive         post_ictal_source_drive
		% return_fields        out_vars

	res.model = parse_model(varargin);
		% dt              electrodes      expansion_rate  grid_size       
		% noise           stim_center     dx              bounds
		% time_constants  SS              IC              K               


%%%% Passed as unmatched parameters from parsed model %%%%
	res.electrodes = parse_electrodes(res.model);
		% centerEC        centerNP        dimsEC          dimsNP
		% scaleEC

	res.time_constants = parse_time(res.model);
		% tau_dD          tau_dVe         tau_dVi         

	res.K = parse_potassium(res.model);
		% k_decay         kD              kR              KtoD
		% KtoVe           KtoVi           tau_K

	res.noise = parse_noise(res.model);
		% noise           noise_sc        noise_sf        


%%%% IC (initial conditions) and SS (steady states) %%%%

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
		
	res.IC = parse_IC(res.model, res.SS);  
		% D11             D22             K               
		% F_ee            F_ei            F_ie            F_ii
		% Phi_ee          Phi_ei          Phi_ie          Phi_ii          
		% phi2_ee         phi2_ei         phi_ee          phi_ei
		% Qe              Qi              Ve              Vi
		% dVe             dVi             map             state

	res.bounds = parse_bounds(res.model);
		% <same as IC>


res.model = check_time_resolution(res.model, res.IC, res.SS);
res.t0 = res.meta.t0_start;

end

function meta = parse_meta(options)
% Save and visualize
[G, p] = get_parser();

p('basename', 'SCM/SCM/SCM');
p('sim_num', []);
p('save', true);  % Save output
p('visualization_rate', 0);  % Show this many frames per second
p('t_step', 1);  % Simulate t_step second intervals
p('t0_start', 0);  % Continue from previous sim (IC will use last from file number t0_start-1)
p('duration', 190);  % Seizure duration
p('padding', [10 30], @(x) isnumeric(x) & numel(x) == 2);  % Padding before and after seizure
p('source_drive', 3);
p('post_ictal_source_drive', 1.5);
p('subsample', Inf);  % Allow subsampling when making mea
p('return_fields', {'Qe', 'Ve'});  % Qe required to make mea
p('out_vars', {'Qe', 'Ve'});  % Define which variables you would like to visualize (can be any in IC)

if isstruct(options), G.parse(options); else, parse(G, options{:}); end
meta = G.Results;

s = 0;
while isempty(meta.sim_num)
	if isempty(dir([meta.basename '_' num2str(s) '*.mat'])), meta.sim_num = s; end
	s = s+1;
end

end

function model = parse_model(options)
% Model parameters

[G, p] = get_parser();

p('stim_center', [0 0], @(x) numel(x) == 2);
p('grid_size', [100 100], @(x) ~all(mod(x, 2)) && numel(x) == 2);  % size of grid to simulate (must be even)
p('dt', 2e-4);
% p('Laplacian', [0 1 0; 1 -4 1; 0 1 0]);  *** switched to local del2
% definition with this in a subfunction
p('dx', .3);  % (cm) (to call each grid point an electrode we need this to be .04);
p('expansion_rate', 1/3, @(x) x >= 0);  % in Hz; set to 0 for fixed source
p('IC', {});  % Initial conditions of the sim (enter options as name-value pairs)
p('SS', {});  % Constants that define the locations of steady states
p('time_constants', {});  % tau parameters controlling rates
p('K', {});  % Potassium (K-related) parameters
p('noise', {});  % Noise parameters
p('electrodes', {});  % electrode positions
p('bounds', {});  % Variables with integration boundaries

% Parse
if isstruct(options), G.parse(options); else, parse(G, options{:}); end
model = G.Results;

% Clean
if isempty(model.time_constants), model.time_constants = G.Unmatched; end
if isempty(model.K), model.K = G.Unmatched; end
if isempty(model.electrodes), model.electrodes = G.Unmatched; end
if isempty(model.noise), model.noise = G.Unmatched; end

% set stimulus center if unset
stim_center = model.stim_center;
if any(stim_center<=0), stim_center = round(model.grid_size * .4); end
model.stim_center = round(stim_center);

end

function bounds = parse_bounds(model)

options = model.bounds;
[G, p] = get_parser(); G.CaseSensitive = true;
validate = @(x) numel(x) == 2 && x(2) > x(1);

p('Dee', [-Inf Inf], validate)  % i <--> i gap-junction diffusive-coupling strength (mm^2)
p('Dii', [.009 1], validate)  % The inhibitory gap junctions cannot pass below a minimum value of 0.009 / dx^2.
p('K',  [-Inf 1], validate)  % extracellular potassium concentration (mm^2)
p('Qe', [-Inf Inf], validate)  % Activity of excitatory population.
p('Qi', [-Inf Inf], validate)  % Activity of inhibitory population.
p('Ve', [-Inf Inf], validate)  % Voltage  of excitatory population.
p('Vi', [-Inf Inf], validate)  % Voltage of inhibitory population.
p('dVe', [-Inf 1.5], validate)  % Excitatory resting potential offset (mV)
p('dVi', [-Inf 0.8], validate)  % Inhibitory resting potential offset (mV)
p('Phi_ee', [-Inf Inf], validate)  % e <--> e synaptic flux
p('Phi_ei', [-Inf Inf], validate)  % e <--> i synaptic flux
p('Phi_ie', [-Inf Inf], validate)  % i <--> e synaptic flux
p('Phi_ii', [-Inf Inf], validate)  % i <--> i synaptic flux
p('phi2_ee', [-Inf Inf], validate)  % Wave dynamics
p('phi2_ei', [-Inf Inf], validate)
p('phi_ee', [-Inf Inf], validate)
p('phi_ei', [-Inf Inf], validate)
p('F_ee', [-Inf Inf], validate)  % flux dynamics
p('F_ei', [-Inf Inf], validate)
p('F_ie', [-Inf Inf], validate)
p('F_ii', [-Inf Inf], validate)

if isstruct(options), G.parse(options), else, G.parse(options{:}); end
bounds = G.Results;
% bounds.D11 = bounds.D1 ./ model.dx.^2;
% bounds.D22 = bounds.D2 ./ model.dx.^2;

end

function noise = parse_noise(model)
%% Noise
[G, p] = get_parser();
options = model.noise;

p('noise_sf', 2);  % noise scale-factor [default: 2, original: 4]
p('noise_sc', 0.2);  % subcortical noise

if isstruct(options), G.parse(options), else, G.parse(options{:}); end
noise = G.Results;

end

function IC = parse_IC(model, SS)
%% Initial conditions
% These are the same as in Waikato-Kramer except that dVe defaults to -1
% instead of 1 (original IC in <map_type='ictal_wavefront'> scenario).
% Additionally, resizing is done via bootstrapping in the case of a grid
% size different from 100x100.

G = inputParser; G.CaseSensitive = true;
p = @(varargin) addParameter(G, varargin{:});  % convenience function

p('map', nan);  % Source of ictal activity (on/off)
p('state', nan);  % Seizure state (ictal/non-ictal)
p('Dee', SS.Dee)  % i <--> i gap-junction diffusive-coupling strength (electrodes)
p('Dii', SS.Dii)  % e <--> e gap-junction diffusive-coupling strength in all space (electrodes)
p('K', 0)  % extracellular potassium concentration (mm^2)
p('Qe', 0)  % Activity of excitatory population.
p('Qi', 0)  % Activity of inhibitory population.
p('Ve', SS.Ve_rest)  % Voltage  of excitatory population.
p('Vi', SS.Vi_rest)  % Voltage of inhibitory population.
p('dVe', -1)  % Excitatory resting potential offset (mV)
p('dVi', .1)  % Inhibitory resting potential offset (mV)
p('Phi_ee', 0)  % e <--> e synaptic flux
p('Phi_ei', 0)  % e <--> i synaptic flux
p('Phi_ie', 0)  % i <--> e synaptic flux
p('Phi_ii', 0)  % i <--> i synaptic flux
p('phi2_ee', 0)  % Wave dynamics
p('phi2_ei', 0)
p('phi_ee', 0)
p('phi_ei', 0)
p('F_ee', 0)  % flux dynamics
p('F_ei', 0)
p('F_ie', 0)
p('F_ii', 0)

if isstruct(model.IC), parse(G, model.IC), else, parse(G, model.IC{:}); end

IC = G.Results;
if isnan(IC.state), [IC.map, IC.state] = update_map(model); end
IC = resize(IC, model);  

end

function model = check_time_resolution(model, IC, SS)
% reduce time resolution if gap-junction diffusive coupling strength (D22)
% or axonal velocity (v) is high

	if model.dt > 2e-4 && (any(IC.Dii(:) >= 0.87) || SS.v(:) > 140)
		model.dt = 2e-4;
		warning('High D2 or HL.v; setting dt to 2e-4.');
	end
	
end

function time_constants = parse_time(model)

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

% Get defaults
SS = SCM_init_globs;  % See this function for details

% ... and assign them to parser
for f = fieldnames(SS)', p(f{:}, SS.(f{:})); end

% Parse
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

% Potassium excitability function
p('k_peak', .5);  % K concentration at which voltage offset is greatest
p('k_width', .1);  % Width of K concentration peak

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
% 
	for f = fieldnames(IC)'
		if numel(IC.(f{:})) == 1
			IC.(f{:}) = IC.(f{:}) * ones(model.grid_size);
		else
			if ~all(size(IC.(f{:})) == model.grid_size)
				error('IC parameter ''%s'' is the wrong size', f{:});
			end
		end
	end

end

function [G, p] = get_parser()
G = inputParser; G.KeepUnmatched = true;
p =@(varargin) addParameter(G, varargin{:});

end