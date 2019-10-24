% default_params
function [res, G] = init_scm_params(varargin)

%% parsing functions
validate = @(x, all) any(validatestring(x, all));  % validate strings
isnonneg = @(x) isnumeric(x) & x >= 0;  % check nonnegative

%% seizing_cortical_field
G = inputParser;
p = @(varargin) addParameter(G, varargin{:});  % convenience function
G.KeepUnmatched = true;

% Save and visualize
p('basename', 'SCM/SCM/SCM');
p('sim_num', 0);
p('SAVE', true);  % Save output
p('visualize_results', false);  %Set this variable to true to create plots during simulation.
p('t_step', 1);  % Simulate t_step second intervals

% Seizure parameters
p('duration', 90);  % Seizure duration
p('padding', [10 30], @(x) isnumeric(x) & numel(x) == 2);  % Padding before and after seizure

% Noise
p('noise', 0.5);             %Noise level

% Model parameters
allMapTypes = {'ictal_wavefront', 'fixed_point_source'};
p('map_type', 'ictal_wavefront', @(x) validate(x, allMapTypes));  % Ictal source 
p('stim_center', [0 0], @(x) numel(x) == 2);
p('IC', {});  % Initial conditions of the sim
p('grid_size', [100 100], @(x) ~mod(x, 2) && numel(x) == 2);  % size of grid to simulate (must be even)


parse(G, varargin{:});
res = G.Results;
if any(res.stim_center<=0), res.stim_center = round(res.grid_size * .4); end
res.stim_center = round(res.stim_center);

%% Initial conditions
G = inputParser;
G.CaseSensitive = true;
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

parse(G, res.IC{:});

res.IC = G.Results;

end