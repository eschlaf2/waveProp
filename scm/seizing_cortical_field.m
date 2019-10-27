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

visualize_results = params.visualize_results;
visualization_rate = params.visualization_rate;
grid_size = params.grid_size;
noise = params.noise;
ictal_wavefront = strcmpi(params.map_type, 'ictal_wavefront');
map = IC.map;
state = IC.state;

%% Parameter
%Parameters for proportion of extracellular potassium.
K0 = 0;         %initial value.
tau_K = 200;    %time-constant (/s).
k_decay = 0.1;  %decay rate (/s).
kD = 1;         %diffusion coefficient (cm^2/s).
KtoVe = 10;     %impact on excitatory population resting voltage.
KtoVi = 10;     %impact on inhibitory population resting voltage.
KtoD  = -50;    %impact on inhibitory gap junction strength.
kR    = 0.15;   %scale reaction term. 
    
tau_dD  = 200;  %inhibitory gap junction time-constant (/s).
tau_dVe = 250;  %excitatory population resting voltage time-constant (/s).
tau_dVi = 250;  %inhibitory population resting voltage time-constant (/s).

% set no. of sampling points (must be even!) along each axis of cortical grid
Nx = grid_size(1);
Ny = grid_size(2);

% dimensions for the cortical grid
dx = spatial_resolution;

% D1 = zeros(Nx,Ny)+0.8/100;              %Set initial i <--> i gap-junction diffusive-coupling strength in all space (cm^2)
% D2 = zeros(Nx,Ny)+0.8;                  %Set initial e <--> e gap-junction diffusive-coupling strength in all space (cm^2)
D2 = IC.D22 * dx.^2;
% K  = zeros(Nx,Ny)+K0;                   %Set initial extracellular ion concentration in all space (cm^2)

% initialize constants for Seizing Cortical Model
global HL
HL = SCM_init_globs;

% define synaptic strengths
rho_e = HL.ge;
rho_i = HL.gi;

% initialize random number generator (from input argument)
rand_state = sum(100*clock);
rng(rand_state, 'v5normal');

noise_sf = 0.2*20*noise;    % noise scale-factor
noise_sc = 0.2;             % subcortical noise

HL.v = 280;                 % axonal conduction velocity (cm/s), [original = 140 cm/s]
HL.Lambda = 4.0;			% inverse-length scale for connectivity (/cm)
HL.gamma_e = 170;           % EPSP decay rate (/s)
HL.gamma_i = 50;            % IPSP decay rate (/s)
HL.tau_e = 0.02;			% excit neuron time-constant (/s) [original = 0.04 s]
HL.tau_i = 0.02;			% inhib neuron time-constant (/s) [original = 0.04 s]

% set time resolution
if all(D2 < 0.87)
	dt = 0.4*1e-3;
else
	dt = 0.2*1e-3;
end

if HL.v > 140
    dt = 0.2*1e-3;
end

% number of time-steps for simulation
Nsteps = round(time_end/dt);
time   = (0:Nsteps-1)'*dt + IC.t0;

% 3x3 Laplacian matrix (used in grid convolution calculations)
Laplacian = [0 1 0; 1 -4 1; 0 1 0];

% diffusion multipliers (these depend on spatial resolution)
% D11 = D1/dx^2;
% D22 = D2/dx^2;

% set up storage vectors and grids
% voltage and activities.

% "Microscale" electrode positions.
xNP = Nx/2;                 
yNP = Ny/2;
halfWin = 4;
win = 2 * halfWin + 1;

% "Macroscale" electrode positions.
rEC = Nx/2 + [      -12, -12, -12,   0,  0,  0,  12, 12, 12]/3;
cEC = Ny/2 + [      -12,   0,  12, -12,  0, 12, -12,  0, 12]/3;

% Output variables for microscale, (*) denotes returned.
 QeNP = zeros(Nsteps,win,win);              %Activity of excitatory population.   (*)
 VeNP = zeros(Nsteps,win,win);              %Voltage  of excitatory population.   (*)
 QiNP = zeros(Nsteps,win,win);              %Activity of inhibitory population.
 ViNP = zeros(Nsteps,win,win);              %Voltage  of inhibitory population.
  DNP = zeros(Nsteps,win,win);              %Inhibitory-to-inhibitory gap junction strength.
  KNP = zeros(Nsteps,win,win);              %Extracellular potassium.
 dVeNP= zeros(Nsteps,win,win);              %Change in resting voltage of excitatory population.
 dViNP= zeros(Nsteps,win,win);              %Change in resting voltage of inhibitory population.

% Output variables for macroscale, (*) denotes returned.
 QeEC = zeros(Nsteps,length(rEC));      %Activity of excitatory population.   (*)
 VeEC = zeros(Nsteps,length(rEC));      %Voltage  of excitatory population.   (*)
 QiEC = zeros(Nsteps,length(rEC));      %Activity of inhibitory population.
 ViEC = zeros(Nsteps,length(rEC));      %Voltage  of inhibitory population.
  DEC = zeros(Nsteps,length(rEC));      %Inhibitory-to-inhibitory gap junction strength.
  KEC = zeros(Nsteps,length(rEC));      %Extracellular potassium.
 dVeEC= zeros(Nsteps,length(rEC));      %Change in resting voltage of excitatory population.
 dViEC= zeros(Nsteps,length(rEC));      %Change in resting voltage of inhibitory population.

%Use as initial conditions the "last" values of previous simulation.
Qe_grid = IC.Qe;
Qi_grid = IC.Qi;

Ve_grid = IC.Ve;
Vi_grid = IC.Vi;

phi_ee = IC.phi_ee;
phi_ei = IC.phi_ei;

phi2_ee = IC.phi2_ee;
phi2_ei = IC.phi2_ei;

F_ee = IC.F_ee;
F_ei = IC.F_ei;
F_ie = IC.F_ie;
F_ii = IC.F_ii;

Phi_ee = IC.Phi_ee;
Phi_ei = IC.Phi_ei;
Phi_ie = IC.Phi_ie;
Phi_ii = IC.Phi_ii;

D11 = IC.D11;
D22 = IC.D22;

del_VeRest = IC.dVe;
del_ViRest = IC.dVi;

K = IC.K;

% noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
B_ee = noise_sf * sqrt(noise_sc* HL.phi_ee_sc / dt);
B_ei = noise_sf * sqrt(noise_sc* HL.phi_ei_sc / dt);

%% Visualization initialization
% Visualization
if visualize_results
	titles = {'Qe', 'Qi', 'K', 'Qe + Qi'};
	clims = {[0 30], [0 30], [0 1], [0 30]};
	tt = -1;
	if ~exist('fig', 'var') || isempty(fig)
		fig = figure;
		h = gobjects(4, 1);
		th = gobjects(4, 1);
		ah = gobjects(4, 1);
		% Image of excitatory population activity.
		  for ii = 1:4
			  h(ii) = subplot(2,2,ii);
			  ah(ii) = imagesc(1:Nx, 1:Ny, zeros(Nx, Ny), clims{ii});
			  th(ii) = title(titles{ii});
			  colormap jet; axis equal; axis tight; axis ij;
		  end

		  % Indicate electrode positions.
		  hold(h(1), 'on')
		  plot(h(1), xNP-1:xNP+1, yNP-1, '*k')
		  plot(h(1), xNP-1:xNP+1, yNP,   '*k')
		  plot(h(1), xNP-1:xNP+1, yNP+1, '*k')
		  for ck=1:length(rEC)
			  plot(h(1), [rEC(ck)-2, rEC(ck)+2, rEC(ck)+2, rEC(ck)-2, rEC(ck)-2], ...
				  [cEC(ck)-2, cEC(ck)-2, cEC(ck)+2, cEC(ck)+2, cEC(ck)-2])
		  end
		  hold(h(1), 'off')
	else
		nax = length(fig.Children);
		for ii = 0:nax - 1
			ah(ii+1) = fig.Children(nax - ii).Children(end);
			th(ii+1) = fig.Children(nax - ii).Title;
		end
	end
end

%% Simulation
for i = 1: Nsteps
	
	if ictal_wavefront && i > 1 && diff(floor(t(i-1:i) / rate))
		[map, state] = update_map(params, state);
	end

    %Save the "microscale" dynamics.
    QeNP(i,:,:) =    Qe_grid(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
    VeNP(i,:,:) =    Ve_grid(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
    QiNP(i,:,:) =    Qi_grid(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
    ViNP(i,:,:) =    Vi_grid(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
     DNP(i,:,:) =        D22(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
     KNP(i,:,:) =          K(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
   dVeNP(i,:,:) = del_VeRest(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
   dViNP(i,:,:) = del_ViRest(xNP-halfWin:xNP+halfWin, yNP-halfWin:yNP+halfWin);
    
    %Save the "macroscale" dynamics.
    for ck=1:length(rEC)
        QeEC(i,ck) = mean(mean(   Qe_grid(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
        VeEC(i,ck) = mean(mean(   Ve_grid(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
        QiEC(i,ck) = mean(mean(   Qi_grid(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
        ViEC(i,ck) = mean(mean(   Vi_grid(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
         DEC(i,ck) = mean(mean(       D22(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
         KEC(i,ck) = mean(mean(         K(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
       dVeEC(i,ck) = mean(mean(del_VeRest(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
       dViEC(i,ck) = mean(mean(del_ViRest(rEC(ck)-1:rEC(ck)+1, cEC(ck)-2:cEC(ck)+1)));
    end    
    
% 1. update wave equations
    phi2_ee_1 = zeros(Nx,Ny);
    phi2_ee_1(2:Nx-1,2:Ny-1) = phi2_ee(2:Nx-1,2:Ny-1) + dt*(-2*HL.v.*HL.Lambda.*phi2_ee(2:Nx-1,2:Ny-1) ...
                             - (HL.v.*HL.Lambda).^2.*phi_ee(2:Nx-1,2:Ny-1) ...
                             + (HL.v.*HL.Lambda).^2.*Qe_grid(2:Nx-1,2:Ny-1))...
                             + dt*(HL.v/dx)^2*convolve2(phi_ee, Laplacian, 'valid');
    phi_ee_1 = phi_ee + dt*phi2_ee;

    phi2_ei_1 = zeros(Nx,Ny);
    phi2_ei_1(2:Nx-1,2:Ny-1) = phi2_ei(2:Nx-1,2:Ny-1) + dt*(-2*HL.v.*HL.Lambda.*phi2_ei(2:Nx-1,2:Ny-1) ...
                             - (HL.v.*HL.Lambda).^2.*phi_ei(2:Nx-1,2:Ny-1) ...
                             + (HL.v.*HL.Lambda).^2.*Qe_grid(2:Nx-1,2:Ny-1)) ...
                             + dt*(HL.v/dx)^2*convolve2(phi_ei, Laplacian, 'valid');
    phi_ei_1 = phi_ei + dt*phi2_ei;

% 2. update the 4 synaptic flux equations (include sc noise)

%%%% E-to-E %%%%
    F_ee_1   = F_ee     +dt*HL.gamma_e.^2*(-2/HL.gamma_e*F_ee - Phi_ee ...
                        +HL.Nee_a*phi_ee ...    %long range
                        +HL.Nee_b*Qe_grid ...   %short range
                        +noise_sc*HL.phi_ee_sc ...    %subcortical (tonic)
                        +B_ee*randn(Nx, Ny));   %subcortical (random)
    Phi_ee_1 = Phi_ee + dt*F_ee;

%%%% E-to-I %%%%
    F_ei_1   = F_ei     +dt*HL.gamma_e.^2*(-2/HL.gamma_e*F_ei - Phi_ei ...
                        +HL.Nei_a*phi_ei ...    %long range
                        +HL.Nei_b*Qe_grid ...   %short range
                        +noise_sc*HL.phi_ei_sc ...    %subcortical (tonic)
                        +B_ei*randn(Nx, Ny));   %subcortical (random)
    Phi_ei_1 = Phi_ei + dt*F_ei;

%%%% I-to-E %%%%
    F_ie_1   = F_ie     +dt*HL.gamma_i.^2*(-2/HL.gamma_i*F_ie - Phi_ie ...
                        +HL.Nie_b*Qi_grid);     %short range
    Phi_ie_1 = Phi_ie + dt*F_ie;

%%%% I-to-I %%%%
    F_ii_1   = F_ii     +dt*HL.gamma_i.^2*(-2/HL.gamma_i*F_ii - Phi_ii ...
                        +HL.Nii_b*Qi_grid);     %short range
    Phi_ii_1 = Phi_ii + dt*F_ii;

% 3. update the soma voltages

    Ve_grid_1 = zeros(Nx,Ny);
    Ve_grid_1(2:Nx-1,2:Ny-1) = Ve_grid(2:Nx-1,2:Ny-1) + dt/HL.tau_e*( (HL.Ve_rest - Ve_grid(2:Nx-1,2:Ny-1)) + del_VeRest(2:Nx-1,2:Ny-1) ...
          + rho_e*Psi_ee(Ve_grid(2:Nx-1,2:Ny-1)).*Phi_ee(2:Nx-1,2:Ny-1) ...      %E-to-E
          + rho_i*Psi_ie(Ve_grid(2:Nx-1,2:Ny-1)).*Phi_ie(2:Nx-1,2:Ny-1) ...      %I-to-E
          + D11(2:Nx-1,2:Ny-1).*convolve2(Ve_grid, Laplacian, 'valid'));

    Vi_grid_1 = zeros(Nx,Ny);
    Vi_grid_1(2:Nx-1,2:Ny-1) = Vi_grid(2:Nx-1,2:Ny-1) + dt/HL.tau_i*( (HL.Vi_rest - Vi_grid(2:Nx-1,2:Ny-1)) + del_ViRest(2:Nx-1,2:Ny-1) ...
          + rho_e*Psi_ei(Vi_grid(2:Nx-1,2:Ny-1)).*Phi_ei(2:Nx-1,2:Ny-1) ...      %E-to-I
          + rho_i*Psi_ii(Vi_grid(2:Nx-1,2:Ny-1)).*Phi_ii(2:Nx-1,2:Ny-1) ...      %I-to-I
          + D22(2:Nx-1,2:Ny-1).*convolve2(Vi_grid, Laplacian, 'valid'));

% 4. update the firing rates
    Qe_grid = HL.Qe_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_e).*(Ve_grid - HL.theta_e)))) ...     %The E voltage must be big enough,
            - HL.Qe_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_e).*(Ve_grid - (HL.theta_e+30)))));   %... but not too big.
    Qi_grid = HL.Qi_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_i).*(Vi_grid - HL.theta_i)))) ...     %The I voltage must be big enough,
            - HL.Qi_max *(1./(1+exp(-pi/(sqrt(3)*HL.sigma_i).*(Vi_grid - (HL.theta_i+30)))));   %... but not too big.

% 5. Update extracellular ion.
    K_1 = zeros(Nx,Ny);
    K_1(2:Nx-1,2:Ny-1) =    K(2:Nx-1,2:Ny-1) ...
                            + dt/tau_K*(-k_decay*K(2:Nx-1,2:Ny-1) ...   % decay term.
                            ...                                         % reaction term.
                            + kR* (Qe_grid(2:Nx-1,2:Ny-1) + Qi_grid(2:Nx-1,2:Ny-1))./(1+exp(-((Qe_grid(2:Nx-1,2:Ny-1) + Qi_grid(2:Nx-1,2:Ny-1))-15))) ...
                            + kD* convolve2(K, Laplacian, 'valid'));     % diffusion term.

% 6. Update inhibitory gap junction strength, and resting voltages.               
    D22_1         = D22        + dt/tau_dD *(KtoD *K);
    del_VeRest_1  = del_VeRest + dt/tau_dVe*(KtoVe*K);
    del_ViRest_1  = del_ViRest + dt/tau_dVi*(KtoVi*K);

    %UPDATE the dynamic variables.
    phi2_ee = phi2_ee_1;
    phi_ee  = phi_ee_1;
    phi2_ei = phi2_ei_1;
    phi_ei  = phi_ei_1;
    
    F_ee   = F_ee_1;
    Phi_ee = Phi_ee_1;
    F_ei   = F_ei_1;
    Phi_ei = Phi_ei_1;
    
    F_ie   = F_ie_1;
    Phi_ie = Phi_ie_1;
    F_ii   = F_ii_1;
    Phi_ii = Phi_ii_1;
    
    Ve_grid = Ve_grid_1;
    Vi_grid = Vi_grid_1;
    
    D22 = max(D22_1,0.1);                   %The inhibitory gap junctions cannot pass below a minimum value of 0.1.
    D11 = D22/100;                          %See definition in [Steyn-Ross et al PRX 2013, Table I].

    del_VeRest = min(del_VeRest_1,1.5);     %The excitatory population resting voltage cannot pass above a maximum value of 1.5.
    
%%%%  Set the "source" locations' excitatory population resting voltage
    ind = find(map==1);
    for k=1:length(ind)
        del_VeRest(ind(k)) = source_del_VeRest;
    end

    del_ViRest = min(del_ViRest_1,0.8);     %The inhibitory population resting voltage cannot pass above a maximum value of 0.8.
    K = min(K_1,1);                         %The extracellular ion cannot pass above a maximum value of 1.0.

%%%%  Implement the no flux boundary conditions.  %%%%

%Top boundary.
phi2_ee(1,:) = phi2_ee(2,:);
 phi_ee(1,:) =  phi_ee(2,:);
phi2_ei(1,:) = phi2_ei(2,:);
 phi_ei(1,:) =  phi_ei(2,:);
   F_ee(1,:) =    F_ee(2,:);
 Phi_ee(1,:) =  Phi_ee(2,:);
   F_ei(1,:) =    F_ei(2,:);
 Phi_ei(1,:) =  Phi_ei(2,:);
   F_ie(1,:) =  F_ie(2,:);
 Phi_ie(1,:) =  Phi_ie(2,:);
   F_ii(1,:) =    F_ii(2,:);
 Phi_ii(1,:) =  Phi_ii(2,:);
Ve_grid(1,:) = Ve_grid(2,:);
Vi_grid(1,:) = Vi_grid(2,:);
Qe_grid(1,:) = Qe_grid(2,:);
Qi_grid(1,:) = Qi_grid(2,:);
    D22(1,:) = D22(2,:);
    D11(1,:) = D11(2,:);
del_VeRest(1,:) = del_VeRest(2,:);
del_ViRest(1,:) = del_ViRest(2,:);
      K(1,:) = K(2,:);

%Bottom boundary.
phi2_ee(Nx,:) = phi2_ee(Nx-1,:);
 phi_ee(Nx,:) =  phi_ee(Nx-1,:);
phi2_ei(Nx,:) = phi2_ei(Nx-1,:);
 phi_ei(Nx,:) =  phi_ei(Nx-1,:);
   F_ee(Nx,:) =    F_ee(Nx-1,:);
 Phi_ee(Nx,:) =  Phi_ee(Nx-1,:);
   F_ei(Nx,:) =    F_ei(Nx-1,:);
 Phi_ei(Nx,:) =  Phi_ei(Nx-1,:);
   F_ie(Nx,:) =    F_ie(Nx-1,:);
 Phi_ie(Nx,:) =  Phi_ie(Nx-1,:);
   F_ii(Nx,:) =    F_ii(Nx-1,:);
 Phi_ii(Nx,:) =  Phi_ii(Nx-1,:);
Ve_grid(Nx,:) = Ve_grid(Nx-1,:);
Vi_grid(Nx,:) = Vi_grid(Nx-1,:);
Qe_grid(Nx,:) = Qe_grid(Nx-1,:);
Qi_grid(Nx,:) = Qi_grid(Nx-1,:);
    D22(Nx,:) = D22(Nx-1,:);
    D11(Nx,:) = D11(Nx-1,:);
del_VeRest(Nx,:) = del_VeRest(Nx-1,:);
del_ViRest(Nx,:) = del_ViRest(Nx-1,:);
      K(Nx,:) = K(Nx-1,:);

%Left boundary.
phi2_ee(:,1) = phi2_ee(:,2);
 phi_ee(:,1) =  phi_ee(:,2);
phi2_ei(:,1) = phi2_ei(:,2);
 phi_ei(:,1) =  phi_ei(:,2);
   F_ee(:,1) =    F_ee(:,2);
 Phi_ee(:,1) =  Phi_ee(:,2);
   F_ei(:,1) =    F_ei(:,2);
 Phi_ei(:,1) =  Phi_ei(:,2);
   F_ie(:,1) =    F_ie(:,2);
 Phi_ie(:,1) =  Phi_ie(:,2);
   F_ii(:,1) =    F_ii(:,2);
 Phi_ii(:,1) =  Phi_ii(:,2);
Ve_grid(:,1) = Ve_grid(:,2);
Vi_grid(:,1) = Vi_grid(:,2);
Qe_grid(:,1) = Qe_grid(:,2);
Qi_grid(:,1) = Qi_grid(:,2);
    D22(:,1) = D22(:,2);
    D11(:,1) = D11(:,2);
del_VeRest(:,1) = del_VeRest(:,2);
del_ViRest(:,1) = del_ViRest(:,2);
      K(:,1) = K(:,2);

%Right boundary.
phi2_ee(:,Ny) = phi2_ee(:,Ny-1);
 phi_ee(:,Ny) =  phi_ee(:,Ny-1);
phi2_ei(:,Ny) = phi2_ei(:,Ny-1);
 phi_ei(:,Ny) =  phi_ei(:,Ny-1);
   F_ee(:,Ny) =    F_ee(:,Ny-1);
 Phi_ee(:,Ny) =  Phi_ee(:,Ny-1);
   F_ei(:,Ny) =    F_ei(:,Ny-1);
 Phi_ei(:,Ny) =  Phi_ei(:,Ny-1);
   F_ie(:,Ny) =    F_ie(:,Ny-1);
 Phi_ie(:,Ny) =  Phi_ie(:,Ny-1);
   F_ii(:,Ny) =    F_ii(:,Ny-1);
 Phi_ii(:,Ny) =  Phi_ii(:,Ny-1);
Ve_grid(:,Ny) = Ve_grid(:,Ny-1);
Vi_grid(:,Ny) = Vi_grid(:,Ny-1);
Qe_grid(:,Ny) = Qe_grid(:,Ny-1);
Qi_grid(:,Ny) = Qi_grid(:,Ny-1);
    D22(:,Ny) = D22(:,Ny-1);
    D11(:,Ny) = D11(:,Ny-1);
del_VeRest(:,Ny) = del_VeRest(:,Ny-1);
del_ViRest(:,Ny) = del_ViRest(:,Ny-1);
      K(:,Ny) = K(:,Ny-1);

      % sanity check!
      if any(any(isnan(Qe_grid)))
          error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
      end
      
      %Visualization
      if visualize_results
          stride2 = 1 / visualization_rate;
		  if floor(time(i)/stride2) > tt, 
			  tt = tt + 1;
			  
              % Image of excitatory population activity.
			  set(ah(1), 'cdata', Qe_grid);
              
              % Image of inhibitory population activity.
			  set(ah(2), 'cdata', Qi_grid);
              
              % Image of extracellular ion proportion.
			  set(ah(3), 'cdata', K)
			  set(th(3), 'string', sprintf('K %2f', mean(K(:))))
              
              % Image of inhibitory gap junction strength.
			  set(ah(4), 'cdata', Qe_grid + Qi_grid);  
              drawnow;
          end
      end
end

%%%% Save the "last" values of the simulation. %%%%

last.Qe = Qe_grid;
last.Qi = Qi_grid;

last.Ve = Ve_grid;
last.Vi = Vi_grid;

last.phi_ee = phi_ee;
last.phi_ei = phi_ei;

last.phi2_ee = phi2_ee;
last.phi2_ei = phi2_ei;

last.F_ee = F_ee;
last.F_ei = F_ei;
last.F_ie = F_ie;
last.F_ii = F_ii;

last.Phi_ee = Phi_ee;
last.Phi_ei = Phi_ei;
last.Phi_ie = Phi_ie;
last.Phi_ii = Phi_ii;

last.D11 = D11;
last.D22 = D22;

last.dVe = del_VeRest;
last.dVi = del_ViRest;

last.K = K;

last.map = map;
last.state = state;

%%%% Define the output variables of simulation.

NP = {};
NP.Qe = single(QeNP);
NP.Ve = single(VeNP);
NP.Qi = single(QiNP);
NP.Vi = single(ViNP);

EC = {};
EC.Qe = single(QeEC);
EC.Ve = single(VeEC);
EC.Qi = single(QiEC);
EC.Vi = single(ViEC);

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