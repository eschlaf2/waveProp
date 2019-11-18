function H = SCM_init_globs
% function H = SCM_init_globs
% Initialize global constants that define the locations of the steady
% states.  These values are stored in the struct H

% 10-Feb-2009 Only initialize constants that determine the steady states;
%		don't initialize those that only alter stability (tau, gamma, v)
%
% Copyright (c) 2009-2012 Alistair Steyn-Ross: asr(at)waikato.ac.nz
%   School of Engineering, Waikato University, NZ 
% This code forms part of the SOinC (slow oscillations in cortex) project
%
% Feb 2009; modified July 2012
%
% Modified by M. A. Kramer for the Seizing Cortical Model (SCM)
% May 2016.

H.tau_e = 0.02;			% excit neuron time-constant (/s) [original = 0.04 s]
H.tau_i = 0.02;			% inhib neuron time-constant (/s) [original = 0.04 s]

% voltage limits
[H.Ve_rev,  H.Vi_rev]  = deal(0, -70);          % reversal potential (mV)
[H.Ve_rest, H.Vi_rest] = deal(-64, -64);        % resting potential (mV)

% gain per synapse at resting voltage (millivolt.sec)
[H.rho_e, H.rho_i] = deal(1.00e-3, -1.05e-3);

% E/IPSP rate constants
H.gamma_e = 170;  % EPSP decay rate (/s)
H.gamma_i = 50;  % IPSP decay rate (/s)

% gap-junction diffusive-coupling strength
H.Dii = 1;          % i <--> i steady state (cm^2) [0, 1]
H.Dee = H.Dii / 100;  % e <--> e initial in all space (cm^2)

% sigmoid characteristics
[H.Qe_max,  H.Qi_max]  = deal(30,    60);       % sigmoid maximum (s^-1)
[H.theta_e, H.theta_i] = deal(-58.5, -58.5);	% sigmoid threshold (mV)
[H.sigma_e, H.sigma_i] = deal(3.0,   5.0);		% sigmoid 'width' (mV)

% connectivities: j-->k convention (dimensionless)
[H.Nee_a, H.Nei_a] = deal(2000, 2000);			% cortico-cortical
[H.Nee_b, H.Nei_b] = deal(800, 800);
[H.Nie_b, H.Nii_b] = deal(600, 600);
% [H.Nee_sc,H.Nei_sc]= deal(50, 50);              % subcortical

% Nee and Nie totals for cortico-cortical plus intracortical
H.Nee_ab = H.Nee_a + H.Nee_b;
H.Nei_ab = H.Nei_a + H.Nei_b;

% default subcortical fluxes
% H.phi_ee_sc = H.Nee_sc * H.Qe_max;  % original [300]
% H.phi_ei_sc = H.Nei_sc * H.Qe_max;  % original [300]
H.phi_ee_sc = 150;  % original [300]
H.phi_ei_sc = 150;  % original [300]

% axonal conduction velocity (cm/s), 
H.v = 280;  % [original = 140 cm/s]

% inverse-length scale for connectivity (/cm)
H.Lambda = 4.0;			

% d/dV derivatives of psi_ij weighting functions
H.d_psi_ee = -1/(H.Ve_rev - H.Ve_rest);
H.d_psi_ei = -1/(H.Ve_rev - H.Vi_rest);
H.d_psi_ie = -1/(H.Vi_rev - H.Ve_rest);
H.d_psi_ii = -1/(H.Vi_rev - H.Vi_rest);


return

end
%------------------------------------------------------------------------
