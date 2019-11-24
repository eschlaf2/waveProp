%   This function is used to simulate a neural field on a two-dimensional
%   plane.  It is an extension of the Waikato Cortical Model originally developed in
%       Moira L. Steyn-Ross, D. A. Steyn-Ross, and J. W. Sleigh
%       Interacting Turing-Hopf Instabilities Drive Symmetry-breaking Transitions
%       in a Mean-field Model of the Cortex: A Mechanism for the Slow Oscillation
%       Phys. Rev. X vol 3(2), e021005 (2013)
%   and available here
%       http://www2.phys.waikato.ac.nz/asr/SteynRoss_papers/
%
%CALLS
% . seizing_cortical_field_IC.m
%
%INPUTS
%  source_del_VeRest = extra offest to resting potential at the source location (mV).
%  map               = matrix of 2-dimensional locations of the source.
%  time_end          = total time of simulation (/s).
%  IC                = structure with initial conditions for model variables. 
%                      There are two options:
%                    = {}     when no initial conditions are specificed.
%                    = last   when using the "last" state of a previous
%                             simulation, see OUTPUT "last".
%OUTPUTS
%  NP       = the activity (Qe) and voltage (Ve) of the excitatory
%             populations at the "microscale".
%  EC       = the activity (Qe) and voltage (Ve) of the excitatory
%             populations at the "macroscale".
%  time     = time axis of the simulation.
%  last     = the value at the last time step of all model variables. Use
%             this output as the input "IC" when running a sequential
%             series of simulations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seizing Cortical Model (a modification of the Waikato Cortical Model)
%   Copyright (c) 2016 M. A. Kramer
%   Department of Mathematics and Statistics, Boston University, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
function [NP, EC, time, last, fig] = seizing_cortical_field(~, time_end, IC, fig, params)

%% Preferences and parameters
if ~exist('params', 'var') || isempty(params), params = init_scm_params(); end

PM = params.meta;
PK = params.K;
M = params.model;
PE = params.electrodes;
PT = params.time_constants;
PN = params.noise;
SS = params.SS;
MX = params.bounds;
PS = params.sigmoids;

%% A few convenience variables

Nx = M.grid_size(1);
Ny = M.grid_size(2);
dx = M.dx;
dt = M.dt;

% initialize random number generator (from input argument)
rand_state = sum(100*clock);
rng(rand_state, 'v5normal');

% number of time-steps for simulation
Nsteps = round(time_end / dt);
time = ( 0 : Nsteps - 1 )' * dt + params.t0;

% noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
B_ee = PN.noise_sf .* sqrt(PN.noise_sc .* SS.phi_ee_sc / dt);
B_ei = PN.noise_sf .* sqrt(PN.noise_sc .* SS.phi_ei_sc / dt);

%% Define output variables

% Output variables, (*) denotes returned.
% Some are visualized, some are returned.

out_vars = PM.out_vars;

% Indices to capture from larger grid for NP and EC
indsNP = get_inds(M.grid_size, PE.centerNP, PE.dimsNP, 1);
indsEC = get_inds(M.grid_size, PE.centerEC, PE.dimsEC, PE.scaleEC);

% Initialize electrode output variables (structs NP, EC)
get_electrode_values;

%% Define dynamic variables
%Use as initial conditions the "last" values of previous simulation.

last = IC;  % initialize
new = last;  % initialize


%% Simulation
new = last;
for ii = 1: Nsteps
	
	% Update equations (update <new> values)
	update_wave_equations;
	update_synaptic_flux_eq;
	update_soma_voltages;
	update_firing_rates;
	update_extracellular_ion;
	update_gap_resting;
	
	% Correct out-of-bounds values
	correct_OOB_values;
	
	% Set the "source" locations' excitatory population resting voltage and
	% expand the wavefront
	update_source

	% sanity check!
	if any(any(isnan(last.Qe)))
		error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
	end
	
	% UPDATE the dynamic variables (pass <new> values to <last>).
	last = new;
	get_electrode_values;
      
	visualize;
end

% Return requested variables.
no_return = out_vars(~ismember(out_vars, PM.return_fields));
NP = rmfield(NP, no_return);
EC = rmfield(EC, no_return);

% ------------------------------------------------------------------------
%% Nested dynamics functions

% 1. update wave equations
	function update_wave_equations
		new.phi2_ee = ...
			last.phi2_ee ...
			+ dt * ( ...
				-2 * SS.v .* SS.Lambda .* last.phi2_ee ...
				- (SS.v .* SS.Lambda).^2 .* last.phi_ee ...
				+ (SS.v .* SS.Lambda).^2 .* last.Qe...
			) ...
			+ dt * (SS.v / dx)^2 * del2_(last.phi_ee);
		 
		new.phi2_ei = last.phi2_ei ...
			+ dt * ( ...
				- 2*SS.v .* SS.Lambda .* last.phi2_ei ...
				- (SS.v .* SS.Lambda).^2 .* last.phi_ei ...
				+ (SS.v .* SS.Lambda).^2 .* last.Qe...
			) ...
			+ dt * (SS.v / dx)^2 * del2_(last.phi_ei);
							 
		new.phi_ee = last.phi_ee + dt * last.phi2_ee;
		new.phi_ei = last.phi_ei + dt * last.phi2_ei;
		
	end

% 2. update the 4 synaptic flux equations (include sc noise)
	function update_synaptic_flux_eq

	%%%% E-to-E %%%%
		new.F_ee = last.F_ee ...
			+ dt * SS.gamma_e.^2 * ( ...
				- 2 / SS.gamma_e * last.F_ee ...
				- last.Phi_ee ...
				+ SS.Nee_a * last.phi_ee ...  % long range
				+ SS.Nee_b * last.Qe ...      % short range
				+ PN.noise_sc .* SS.phi_ee_sc ... % subcortical (tonic)
				+ B_ee .* randn(Nx, Ny) ...       % subcortical (random)
			);
		
		new.Phi_ee = last.Phi_ee + dt * last.F_ee;

	%%%% E-to-I %%%%
		new.F_ei = last.F_ei ...
			+ dt * SS.gamma_e.^2 * ( ...
				- 2 / SS.gamma_e * last.F_ei ...
				- last.Phi_ei ...
				+ SS.Nei_a * last.phi_ei ...    %long range
				+ SS.Nei_b * last.Qe ...   %short range
				+ PN.noise_sc .* SS.phi_ei_sc ...    %subcortical (tonic)
				+ B_ei .* randn(Nx, Ny)...   %subcortical (random)
			);
		
		new.Phi_ei = last.Phi_ei + dt * last.F_ei;

	%%%% I-to-E %%%%
		new.F_ie = last.F_ie ...
			+ dt * SS.gamma_i.^2 * ( ...
				- 2 / SS.gamma_i * last.F_ie ...
				- last.Phi_ie ...
				+ SS.Nie_b * last.Qi ...     %short range
			);
		
		new.Phi_ie = last.Phi_ie + dt * last.F_ie;

	%%%% I-to-I %%%%
		new.F_ii = last.F_ii ...
			+ dt * SS.gamma_i.^2 * ( ...
				- 2 / SS.gamma_i * last.F_ii ...
				- last.Phi_ii ...
				+ SS.Nii_b * last.Qi ...  %short range
			);
		new.Phi_ii = last.Phi_ii + dt * last.F_ii;

	end

% 3. update the soma voltages
	function update_soma_voltages
		new.Ve = last.Ve ...
			+ dt / SS.tau_e * ( ...
				(SS.Ve_rest - last.Ve) ...
				+ last.dVe ...
				+ SS.rho_e * Psi_ee(last.Ve) .* last.Phi_ee ...      %E-to-E
                + SS.rho_i * Psi_ie(last.Ve) .* last.Phi_ie ...      %I-to-E
				+ last.Dee/dx^2 .* del2_(last.Ve) ...
			);
		new.Vi = last.Vi ...
			+ dt / SS.tau_i * ( ...
				(SS.Vi_rest - last.Vi) ...
				+ last.dVi ...
				+ SS.rho_e * Psi_ei(last.Vi) .* last.Phi_ei ...      %E-to-I
				+ SS.rho_i * Psi_ii(last.Vi) .* last.Phi_ii ...      %I-to-I
				+ last.Dii/dx^2 .* del2_(last.Vi) ...
			);
	end

% 4. update the firing rates
	function update_firing_rates
		
		new.Qe = ...
			SS.Qe_max * ( ...
				1 ./ ( 1 + ...
					exp(-pi / ( sqrt(3) * SS.sigma_e ) .* ( last.Ve - SS.theta_e )) ...
				) ...
			) ...     %The E voltage must be big enough,
			- SS.Qe_max * ( ...
				1 ./ ( 1 + ...
					exp(-pi / ( sqrt(3) * SS.sigma_e ) .* ( last.Ve - ( SS.theta_e + 30 ) )) ...
				) ...
			);   %... but not too big.
			
		new.Qi = ...
			SS.Qi_max * ( ...
				1 ./ (1 + ...
					exp(-pi / ( sqrt(3) * SS.sigma_i ) .* ( last.Vi - SS.theta_i )) ...
				) ...
			) ...     %The I voltage must be big enough,
			- SS.Qi_max * ( ...
				1 ./ ( 1 + ...
					exp(-pi / ( sqrt(3) * SS.sigma_i ) .* ( last.Vi - ( SS.theta_i + 30 ) )) ...
				) ...
			);   %... but not too big.
		
	end

% 5. Update extracellular ion.
	function update_extracellular_ion
		new.K = ...
			last.K ...
			+ dt / PK.tau_K .* ( ...
				- PK.k_decay .* last.K ...  % decay term.
				+ PK.kR .* ...              % reaction term.
					1 ./ ( 1 + exp( -( last.Qe + last.Qi ) + 15) ) ...
			) ...
			+ dt * PK.kD./dx^2 * del2_(last.K);  % diffusion term.
			
	end

% 6. Update inhibitory gap junction strength, and resting voltages.  
	function update_gap_resting
		new.Dii = last.Dii + dt / PT.tau_dD * ( PK.KtoD .* wD(last.K) - (last.Dii - SS.Dii));
		new.Dee = last.Dii/100;                %See definition in [Steyn-Ross et al PRX 2013, Table I].
		new.dVe = last.dVe + dt / PT.tau_dVe .* ( PK.KtoVe .* wdVe(last.K) - last.dVe);
		new.dVi = last.dVi + dt / PT.tau_dVi .* ( PK.KtoVi .* wdVe(last.K) - last.dVi);
	end

% Correct out of bounds values
	function correct_OOB_values
		
		for f = fieldnames(last)'
			f = f{:}; %#ok<FXSET>
			if ismember(f, {'map', 'state'}) || all(isinf(MX.(f))), continue, end
			new.(f) = min(max(new.(f), MX.(f)(1)), MX.(f)(2));
		end
		
	end

% Update source and expand wavefront
	function update_source
		if time(ii) > 0 
			if time(ii) <= PM.duration
				[new.map, new.state] = update_map(last.state, M.expansion_rate * dt / dx^2, M.excitability_map);
				new.dVe(new.map) = PM.source_drive; 
				if isempty(PM.source), return, end
				which_source = mod(floor(time(ii) / 2), size(PM.source, 3)) + 1;
				new.dVe(PM.source(:, :, which_source)) = PM.source_drive;
			else
				new.map = last.map; new.state = last.state;
				if isnan(PM.post_ictal_source_drive), return; end			
				new.dVe(new.map) = PM.post_ictal_source_drive;
			end
		end
	end

%% Nested logistical functions

% Store electrode values
	function get_electrode_values
		
		if ~exist('NP', 'var')  % Initialize
			for v = out_vars
				NP.(v{:}) = zeros([Nsteps, PE.dimsNP], 'single');  % Microscale
				EC.(v{:}) = zeros([Nsteps, PE.dimsEC], 'single');  % Macroscale
			end
		else  % ... or store
			for v = out_vars 
				NP.(v{:})(ii, :, :) = reshape(last.(v{:})(indsNP), [1, PE.dimsNP]); 
				
				% Take the local mean for ECOG electrodes
				blurEC = conv2(last.(v{:}), ones(3) ./ 9, 'same');  
				EC.(v{:})(ii, :, :) = reshape(blurEC(indsEC), [1, PE.dimsEC]);
			end
		end
	end

% Visualize results
	function visualize
		if PM.visualization_rate > 0
			if ~exist('fig', 'var') || isempty(fig)
				fig = create_fig(out_vars, M.grid_size, indsNP, indsEC);
			end
			if diff(floor((time(ii) - [dt 0]) * PM.visualization_rate))
				for sp = 1:numel(fig.ih)
					set(fig.ih(sp), 'cdata', last.(out_vars{sp}))
					colorbar
				end

				set(fig.ah, 'string', sprintf('T = %0.3f', time(ii)));
				drawnow;
			end
		end
	end

%% Nested weighting functions

% Excitability (voltage offset) as a function of potassium
	function y = wdVe(K)
		y = 2 ./ ( 1 + exp( -( 2/PS.kdVe_width * (K - PS.kdVe_center) ) ) ) - 1;
	end
	
% Gapjunction functionality as a function of potassium
	function y = wD(K)
		y = 1 ./ (1 + exp( -( 2/PS.kD_width * (K - PS.kD_center) ) ));
	end

% e-to-e reversal-potential weighting function
	function weight = Psi_ee(V)
		weight = (SS.Ve_rev - V)./(SS.Ve_rev - SS.Ve_rest);
	end

% e-to-i reversal-potential weighting function
	function weight = Psi_ei(V)
		weight = (SS.Ve_rev - V)./(SS.Ve_rev - SS.Vi_rest);
	end

% i-to-e reversal-potential weighting function
	function weight = Psi_ie(V)
		weight = (SS.Vi_rev - V)./(SS.Vi_rev - SS.Ve_rest);
	end

% i-to-i reversal potential weighting function
	function weight = Psi_ii(V)
		weight = (SS.Vi_rev - V)./(SS.Vi_rev - SS.Vi_rest);
	end

end

% ------------------------------------------------------------------------
%% Supplementary functions

%------------------------------------------------------------------------
function Y = del2_(X)

% L = [0 1 0; 1 -4 1; 0 1 0];  % 5-point stencil Laplacian
% r = 4/3 * (1 + 1 / sqrt(2));
% b = 1 / r; a = 1 / (sqrt(2) * r);
a = .25; b = .5;
L = [a b a; b -3 b; a b a];  % 9-point stencil Laplacian

% zero-flux BCs
	X = [repmat(X(1, :), 2, 1); X; repmat(X(end, :), 2, 1)];
	X = [repmat(X(:, 1), 1, 2), X, repmat(X(:, end), 1, 2)];

	Y = conv2(X, L, 'valid');
	Y = Y(2:end-1, 2:end-1);
end

%------------------------------------------------------------------------
function fig = create_fig(out_vars, grid_size, addyNP, addyEC)
	stop_at = 'seizing_cortical_field>update_gap_resting';
	titles = out_vars;
	N = numel(titles);
% 	clims = {[0 30], [0 30], [0 1], [0 30]};
	fig.fig = figure;
	fig.quit_early = false;
	% Create 'Quit' pushbutton in figure window
	uicontrol('units','normal','position',[.45 .02 .13 .07], ...
	    'callback',['dbstop in ' stop_at],...
	    'fontsize',10,'string','Debug');
	dbclear seizing_cortical_field>visualize
	h = gobjects(N, 1);
	th = gobjects(N, 1);
	ih = gobjects(N, 1);
	% Image of excitatory population activity.
	nr = floor(sqrt(N));
	nc = ceil(N / nr);
	for ii = 1:N
		h(ii) = subplot(nr,nc,ii);
		ih(ii) = imagesc(zeros(grid_size));
		th(ii) = title(titles{ii});
		colormap bone; axis equal; axis tight; axis ij; colorbar;
	end

	% Indicate electrode positions.
	colorEC = [0 1 1];
	colorNP = [1 0 1];
	
	hold(h(1), 'on')
	
	[xE, yE] = ind2sub(grid_size, addyEC);
	[xN, yN] = ind2sub(grid_size, addyNP);
	kE = convhull(xE(:), yE(:));
	kN = convhull(xN(:), yN(:));
	plot(h(1), xE(kE), yE(kE), '-', 'color', colorEC);
	plot(h(1), xN(kN), yN(kN), '-', 'color', colorNP);
	legend(h(1), {'EC', 'NP'}, 'location', 'southeast')
	legend(h(1), 'boxoff')

	hold(h(1), 'off')
	ah = annotation('textbox', [0 .99 0 0], 'string', 'T = ', ...
		'FitBoxToText', 'on', 'LineStyle', 'none');
	
	fig.h = h;
	fig.ih = ih;
	fig.th = th;
	fig.ah = ah;
end

%------------------------------------------------------------------------
function inds = get_inds(grid_size, center, dims, scale)

	x_offset = ( (1:dims(1)) - floor(dims(1)/2) ) * scale;
	y_offset = ( (1:dims(2)) - floor(dims(2)/2) ) * scale;
	[xx, yy] = ndgrid(x_offset, y_offset);

	inds = sub2ind(grid_size, ...
		center(1) + xx(:), ...
		center(2) + yy(:));
end

%------------------------------------------------------------------------

