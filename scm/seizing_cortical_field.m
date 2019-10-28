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
function [NP, EC, time, last, fig] = seizing_cortical_field(source_del_VeRest, time_end, IC, fig, params)

%% Preferences and parameters
if ~exist('params', 'var') || isempty(params), params = init_scm_params(), end

PM = params.meta;
PK = params.potassium;
IC = params.IC;
M = params.model;
PE = params.electrodes;
PT = params.time_constants;
PN = params.noise;

% global HL
HL = params.HL;

% struct2var(params.meta);
% struct2var(params.potassium);

%% A few convenience variables

ictal_wavefront = strcmpi(M.map_type, 'ictal_wavefront');
Nx = M.grid_size(1);
Ny = M.grid_size(2);
dx = M.spatial_resolution;
dt = M.dt;

% initialize random number generator (from input argument)
rand_state = sum(100*clock);
rng(rand_state, 'v5normal');

% number of time-steps for simulation
Nsteps = round(time_end/dt);
time   = (0:Nsteps-1)'*dt + IC.t0;

%% Define output variables
% Output variables, (*) denotes returned.
out_vars = {...
	'Qe', ...  %Activity of excitatory population.   (*)
	'Ve', ...  %Voltage  of excitatory population.   (*)
	'Qi', ...  %Activity of inhibitory population.
	'Vi', ...  %Voltage of inhibitory population.
	'D', ...   %Inhibitory-to-inhibitory gap junction strength.
	'K', ...   %Extracellular potassium.
	'dVe', ... %Change in resting voltage of excitatory population.
	'dVi', ... %Change in resting voltage of inihibitory population.
	};

addyNP = get_inds(M.grid_size, PE.centerNP, PE.dimsNP, 1);
addyEC = get_inds(M.grid_size, PE.centerEC, PE.dimsEC, PE.scaleEC);

% Output variables for microscale
for v = out_vars
	NP.(v{:}) = zeros(Nsteps, PE.dimsNP(1), PE.dimsNP(2), 'single'); 
end

% Output variables for macroscale
for v = out_vars
	EC.(v{:}) = zeros(Nsteps, PE.dimsEC(1), PE.dimsEC(2), 'single'); 
end

%% Define dynamic variables
%Use as initial conditions the "last" values of previous simulation.
dynamic_vars = fieldnames(IC)';

% noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
B_ee = PN.noise_sf * sqrt(PN.noise_sc* HL.phi_ee_sc / dt);
B_ei = PN.noise_sf * sqrt(PN.noise_sc* HL.phi_ei_sc / dt);

%% Visualization initialization
% Visualization
if PM.visualize_results, visualization(fig, M.grid_size, addyNP, addyEC); end

%% Simulation
last = IC;
for i = 1: Nsteps
	
	if ictal_wavefront && diff(floor(time(i) - [dt 0] / M.expansion_rate))
		[last.map, last.state] = update_fun(last.state);
	end

    %Save the "microscale" dynamics.
	for v = out_vars 
		NP.(v{:})(i, :, :) = last.(v{:})(addyNP); 
	end
    
    %Save the "macroscale" dynamics.
	for v = out_vars 
		blurEC = conv2(last.(v{:}), ones(3), 'same') ./ 9;
		EC.(v{:})(i, :, :) = blurEC;
	end
    
% 1. update wave equations
    new.phi2_ee = zeros(Nx,Ny);
    new.phi2_ee(2:Nx-1,2:Ny-1) = last.phi2_ee(2:Nx-1,2:Ny-1) + dt*(-2*HL.v.*HL.Lambda.*last.phi2_ee(2:Nx-1,2:Ny-1) ...
                             - (HL.v.*HL.Lambda).^2.*last.phi_ee(2:Nx-1,2:Ny-1) ...
                             + (HL.v.*HL.Lambda).^2.*last.Qe(2:Nx-1,2:Ny-1))...
                             + dt*(HL.v/dx)^2*convolve2(last.phi_ee, M.Laplacian, 'valid');
    new.phi_ee = last.phi_ee + dt*last.phi2_ee;

    new.phi2_ei = zeros(Nx,Ny);
    new.phi2_ei(2:Nx-1,2:Ny-1) = last.phi2_ei(2:Nx-1,2:Ny-1) + dt*(-2*HL.v.*HL.Lambda.*last.phi2_ei(2:Nx-1,2:Ny-1) ...
                             - (HL.v.*HL.Lambda).^2.*last.phi_ei(2:Nx-1,2:Ny-1) ...
                             + (HL.v.*HL.Lambda).^2.*last.Qe(2:Nx-1,2:Ny-1)) ...
                             + dt*(HL.v/dx)^2*convolve2(last.phi_ei, M.Laplacian, 'valid');
    new.phi_ei = last.phi_ei + dt*last.phi2_ei;

% 2. update the 4 synaptic flux equations (include sc noise)

%%%% E-to-E %%%%
    new.F_ee   = last.F_ee +dt*HL.gamma_e.^2*(-2/HL.gamma_e*last.F_ee - last.Phi_ee ...
                        +HL.Nee_a*last.phi_ee ...    %long range
                        +HL.Nee_b*last.Qe ...   %short range
                        +noise_sc*HL.phi_ee_sc ...    %subcortical (tonic)
                        +B_ee*randn(Nx, Ny));   %subcortical (random)
    new.Phi_ee = last.Phi_ee + dt*last.F_ee;

%%%% E-to-I %%%%
    new.F_ei   = last.F_ei+dt*HL.gamma_e.^2*(-2/HL.gamma_e*last.F_ei - last.Phi_ei ...
                        +HL.Nei_a*last.phi_ei ...    %long range
                        +HL.Nei_b*last.Qe ...   %short range
                        +PN.noise_sc*HL.phi_ei_sc ...    %subcortical (tonic)
                        +B_ei*randn(Nx, Ny));   %subcortical (random)
    new.Phi_ei = last.Phi_ei + dt*last.F_ei;

%%%% I-to-E %%%%
    new.F_ie   = last.F_ie     +dt*HL.gamma_i.^2*(-2/HL.gamma_i*last.F_ie - last.Phi_ie ...
                        +HL.Nie_b*last.Qi);     %short range
    new.Phi_ie = last.Phi_ie + dt*last.F_ie;

%%%% I-to-I %%%%
    new.F_ii   = last.F_ii     +dt*HL.gamma_i.^2*(-2/HL.gamma_i*last.F_ii - last.Phi_ii ...
                        +HL.Nii_b*last.Qi);     %short range
    new.Phi_ii = last.Phi_ii + dt*last.F_ii;

% 3. update the soma voltages

    new.Ve_grid = zeros(Nx,Ny);
    new.Ve_grid(2:Nx-1,2:Ny-1) = last.Ve(2:Nx-1,2:Ny-1) + dt/HL.tau_e*( (HL.Ve_rest - last.Ve(2:Nx-1,2:Ny-1)) + last.dVe(2:Nx-1,2:Ny-1) ...
          + HL.rho_e*Psi_ee(last.Ve(2:Nx-1,2:Ny-1)).*last.Phi_ee(2:Nx-1,2:Ny-1) ...      %E-to-E
          + HL.rho_i*Psi_ie(last.Ve(2:Nx-1,2:Ny-1)).*last.Phi_ie(2:Nx-1,2:Ny-1) ...      %I-to-E
          + last.D11(2:Nx-1,2:Ny-1).*convolve2(last.Ve, M.Laplacian, 'valid'));

    new.Vi_grid = zeros(Nx,Ny);
    new.Vi_grid(2:Nx-1,2:Ny-1) = last.Vi(2:Nx-1,2:Ny-1) + dt/HL.tau_i*( (HL.Vi_rest - last.Vi(2:Nx-1,2:Ny-1)) + last.dVi(2:Nx-1,2:Ny-1) ...
          + HL.rho_e*Psi_ei(last.Vi(2:Nx-1,2:Ny-1)).*last.Phi_ei(2:Nx-1,2:Ny-1) ...      %E-to-I
          + HL.rho_i*Psi_ii(last.Vi(2:Nx-1,2:Ny-1)).*last.Phi_ii(2:Nx-1,2:Ny-1) ...      %I-to-I
          + last.D22(2:Nx-1,2:Ny-1).*convolve2(last.Vi, M.Laplacian, 'valid'));

% 4. update the firing rates
    new.Qe = HL.Qe_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_e).*(last.Ve - HL.theta_e)))) ...     %The E voltage must be big enough,
            - HL.Qe_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_e).*(last.Ve - (HL.theta_e+30)))));   %... but not too big.
    new.Qi = HL.Qi_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_i).*(last.Vi - HL.theta_i)))) ...     %The I voltage must be big enough,
            - HL.Qi_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_i).*(last.Vi - (HL.theta_i+30)))));   %... but not too big.

% 5. Update extracellular ion.
    new.K = zeros(Nx,Ny);
    new.K(2:Nx-1,2:Ny-1) =    last.K(2:Nx-1,2:Ny-1) ...
                            + dt/PT.tau_K*(-PK.k_decay*last.K(2:Nx-1,2:Ny-1) ...   % decay term.
                            ...                                         % reaction term.
                            + PK.kR* (new.Qe(2:Nx-1,2:Ny-1) + new.Qi(2:Nx-1,2:Ny-1))./(1+exp(-((new.Qe(2:Nx-1,2:Ny-1) + new.Qi(2:Nx-1,2:Ny-1))-15))) ...
                            + PK.kD* convolve2(last.K, M.Laplacian, 'valid'));     % diffusion term.

% 6. Update inhibitory gap junction strength, and resting voltages.               
    new.D22         = last.D22        + dt/PT.tau_dD *(PK.KtoD *last.K);
    new.dVe  = last.dVe + dt/PT.tau_dVe*(PK.KtoVe*last.K);
    new.dVi  = last.dVi + dt/PT.tau_dVi*(PK.KtoVi*last.K);

    %UPDATE the dynamic variables.
	for v = dynamic_vars
		last.(v{:}) = new.(v{:});
	end
    
	% Correct out of bounds values
    last.D22 = max(last.D22,0.1);                   %The inhibitory gap junctions cannot pass below a minimum value of 0.1.
    last.D11 = last.D22/100;                          %See definition in [Steyn-Ross et al PRX 2013, Table I].

    last.dVe = min(last.dVe,1.5);     %The excitatory population resting voltage cannot pass above a maximum value of 1.5.
	last.dVi = min(last.dVi,0.8);     %The inhibitory population resting voltage cannot pass above a maximum value of 0.8.
	last.K = min(last.K,1);                         %The extracellular ion cannot pass above a maximum value of 1.0.

%%%%  Set the "source" locations' excitatory population resting voltage
	last.dVe(last.map) = source_del_VeRest;

%%%%  Implement the no flux boundary conditions.  %%%%
	for v = dynamic_vars
		temp = last.(v{:});
		last.(v{:}) = padarray(temp(2:end-1, 2:end-1), [1 1], 'replicate');
	end

	% sanity check!
	if any(any(isnan(last.Qe)))
		error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
	end
      
	%Visualization
	if PM.visualize_results && floor(time(i) * PM.visualization_rate) > tt 
		visualization_update(ah, th, Qe, Qi, K),
		tt = tt + 1;
	end
end

% if ictal_wavefront, last.state = state; end

%%%% Define the output variables of simulation.

NP = rmfield(NP, {'D', 'K', 'dVe', 'dVi'});
EC = rmfield(EC, {'D', 'K', 'dVe', 'dVi'});

end

%------------------------------------------------------------------------
function visualization_update(ah, th, Qe, Qi, K)
% Image of excitatory population activity.
set(ah(1), 'cdata', Qe);

% Image of inhibitory population activity.
set(ah(2), 'cdata', Qi);

% Image of extracellular ion proportion.
set(ah(3), 'cdata', K)
set(th(3), 'string', sprintf('K %2f', mean(K(:))))

% Image of inhibitory gap junction strength.
set(ah(4), 'cdata', Qe + Qi);  
drawnow;
end

%------------------------------------------------------------------------
function [fig, ah, th] = visualization(fig, grid_size, addyNP, addyEC)
	titles = {'Qe', 'Qi', 'K', 'Qe + Qi'};
	clims = {[0 30], [0 30], [0 1], [0 30]};
	if ~exist('fig', 'var') || isempty(fig)
		fig = figure;
		h = gobjects(4, 1);
		th = gobjects(4, 1);
		ah = gobjects(4, 1);
		% Image of excitatory population activity.
		  for ii = 1:4
			  h(ii) = subplot(2,2,ii);
			  ah(ii) = imagesc(zeros(grid_size), clims{ii});
			  th(ii) = title(titles{ii});
			  colormap jet; axis equal; axis tight; axis ij;
		  end

		  % Indicate electrode positions.
		  hold(h(1), 'on')
		  map = nan(grid_size);
		  map(addyNP) = 1;
		  map(addyEC) = 2;
		  h = pcolor(map);
		  h.FaceAlpha = .5;
		  h.LineStyle = 'none';
		  
		  hold(h(1), 'off')
	else
		nax = length(fig.Children);
		for ii = 0:nax - 1
			ah(ii+1) = fig.Children(nax - ii).Children(end);
			th(ii+1) = fig.Children(nax - ii).Title;
		end
	end
end

%------------------------------------------------------------------------
function inds = get_inds(grid_size, center, dims, scale)

x_offset = (1:dims(1) - floor(dims(1) / 2)) * scale;
y_offset = (1:dims(2) - floor(dims(2) / 2)) * scale;

inds = sub2ind(grid_size, ...
	center(1) + x_offset, ...
	center(2) + y_offset);

end

%------------------------------------------------------------------------
function weight = Psi_ee(V)
% e-to-e reversal-potential weighting function

global HL
weight = (HL.Ve_rev - V)/(HL.Ve_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ei(V)
% e-to-i reversal-potential weighting function

global HL
weight = (HL.Ve_rev - V)/(HL.Ve_rev - HL.Vi_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ie(V)
% i-to-e reversal-potential weighting function

global HL
weight = (HL.Vi_rev - V)/(HL.Vi_rev - HL.Ve_rest);
end

%------------------------------------------------------------------------
function weight = Psi_ii(V)
% i-to-i reversal potential weighting function

global HL
weight = (HL.Vi_rev - V)/(HL.Vi_rev - HL.Vi_rest);
end

%------------------------------------------------------------------------