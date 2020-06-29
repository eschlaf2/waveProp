classdef SCM < handle
% classdef SCM
% Initialize global constants that define the locations of the steady
% states.  

% 10-Feb-2009 Only initialize constants that determine the steady states
%		don't initialize those that only alter stability (tau, gamma, v)
%
% Copyright (c) 2009-2012 Alistair Steyn-Ross: asr(at)waikato.ac.nz
%   School of Engineering, Waikato University, NZ 
% This code forms part of the SOinC (slow oscillations in cortex) project
%
% Feb 2009 modified July 2012
%
% Modified by M. A. Kramer for the Seizing Cortical Model (SCM)
% May 2016.


%% Methods

	methods
		function S = SCM(obj, grid_size)
			if nargin < 1, return; end
			if nargin < 2, grid_size = []; end
			
			if ~isempty(obj)
				for ff = string(fieldnames(obj)'), S.(ff) = obj.(ff); end
			end
			if ~isempty(grid_size), S.resize(grid_size); end
			
		end

	end

%% Current state (IC)
	properties 
		%% Initial conditions
			% These are the same as in Waikato-Kramer except that dVe defaults to -1
			% instead of 1 (original IC in <map_type='ictal_wavefront'> scenario).

			
			Dee = 0.01;  % e <--> e gap-junction diffusive-coupling strength (electrodes)
			Dii = 1;  % i <--> i gap-junction diffusive-coupling strength in all space (electrodes)
			K = 0;  % extracellular potassium concentration (cm^2)
			Qe = 0;  % Activity of excitatory population.
			Qi = 0;  % Activity of inhibitory population.
			Ve = -64 % Voltage  of excitatory population.
			Vi = -64 % Voltage of inhibitory population.
			dVe = -1;  % Excitatory resting potential offset (mV)
			dVi = .1;  % Inhibitory resting potential offset (mV)
			Phi_ee = 0;  % e <--> e synaptic flux
			Phi_ei = 0;  % e <--> i synaptic flux
			Phi_ie = 0;  % i <--> e synaptic flux
			Phi_ii = 0;  % i <--> i synaptic flux
			phi2_ee = 0;  % Wave dynamics
			phi2_ei = 0;
			phi_ee = 0;
			phi_ei = 0;
			F_ee = 0;  % flux dynamics
			F_ei = 0;
			F_ie = 0;
			F_ii = 0;
			
			% Seizure source and state
			map logical
			state double
					

	end
	
	methods
		
		function resize(S, grid_size)
		% resize(SCM, grid_size)
		% Converts scalar to matrix
		
% 			fields = metaclass(S).PropertyList;
% 			mask = ~cat(1, fields.Constant) & ~cat(1, fields.Dependent);
% 			for ff = string({fields(mask).Name})
			for ff = string(properties(S)')
				if ismember(ff, {'map', 'state'}), continue; end
				if all(size(S.(ff)) == grid_size), continue; end
				S.(ff) = S.(ff)(1) .* ones(grid_size); 
			end
		end
		
	end
end


%%
%% Steady states
% 	properties  (Constant = true)
% 		tau_e = 0.02			% excit neuron time-constant (/s) [original = 0.04 s]
% 		tau_i = 0.02			% inhib neuron time-constant (/s) [original = 0.04 s]
% 
% 		% voltage limits
% 		Ve_rev = 0          % reversal potential (mV)
% 		Vi_rev = -70
% 		Ve_rest = -64        % resting potential (mV)
% 		Vi_rest = -64
% 
% 		% gain per synapse at resting voltage (millivolt.sec)
% 		rho_e = 1.00e-3
% 		rho_i = -1.05e-3
% 
% 		% E/IPSP rate constants
% 		gamma_e = 170  % EPSP decay rate (/s)
% 		gamma_i = 50  % IPSP decay rate (/s)
% 
% 		% gap-junction diffusive-coupling strength
% 		D = 1          % i <--> i steady state (cm^2) [0, 1]
% 
% 		% sigmoid characteristics
% 		Qe_max = 30       % sigmoid maximum (s^-1)
% 		Qi_max = 60
% 		theta_e = -58.5	 % sigmoid threshold (mV)
% 		theta_i = -58.5
% 		sigma_e = 3.0	% sigmoid 'width' (mV)
% 		sigmi_i = 5.0
% 
% 		% connectivities: j-->k convention (dimensionless)			
% 		Nee_a = 2000  % cortico-cortical
% 		Nei_a = 2000
% 		Nee_b = 800
% 		Nei_b = 800
% 		Nie_b = 600
% 		Nii_b = 600
% 		% [Nee_sc,Nei_sc]= deal(50, 50)              % subcortical
% 
% 		
% 		% default subcortical fluxes
% 		% phi_ee_sc = Nee_sc * Qe_max  % original [300]
% 		% phi_ei_sc = Nei_sc * Qe_max  % original [300]
% 		phi_ee_sc = 150  % original [300]
% 		phi_ei_sc = 150  % original [300]
% 
% 		% axonal conduction velocity (cm/s), 
% 		v = 280  % [original = 140 cm/s]
% 
% 		% inverse-length scale for connectivity (/cm)
% 		Lambda = 4.0			
% 	end
% 	
% 	properties (Dependent = true)
% 		% d/dV derivatives of psi_ij weighting functions
% 		d_psi_ee 
% 		d_psi_ei 
% 		d_psi_ie 
% 		d_psi_ii 
% 		
% 		% Nee and Nie totals for cortico-cortical plus intracortical
% 		Nee_ab
% 		Nei_ab
% 
% 	end
% 	
% 	methods
% 		function d_psi_ee = get.d_psi_ee(S)
% 			d_psi_ee = -1/(S.Ve_rev - S.Ve_rest);
% 		end
% 		function d_psi_ei = get.d_psi_ei(S)
% 			d_psi_ei = -1/(S.Ve_rev - S.Vi_rest);
% 		end
% 		function d_psi_ie = get.d_psi_ie(S)
% 			d_psi_ie = -1/(S.Vi_rev - S.Ve_rest);
% 		end
% 		function d_psi_ii = get.d_psi_ii(S)
% 			d_psi_ii = -1/(S.Vi_rev - S.Vi_rest);
% 		end
% 		
% 		% Nee and Nie totals for cortico-cortical plus intracortical
% 		function N = get.Nee_ab(S)
% 			N = S.Nee_a + S.Nee_b;
% 		end
% 		function N = get.Nei_ab(S)
% 			N = S.Nei_a + S.Nei_b;
% 		end
% 
% 	end
	


