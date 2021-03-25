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
function [NP, EC, time, last, fig] = seizing_cortical_field(~, time_end, IC, fig, scm)

%% Preferences and parameters
if ~exist('scm', 'var') || isempty(scm), scm = SCM; end

%% A few convenience variables

Nx = scm.grid_size(1);
Ny = scm.grid_size(2);
dx = scm.dx;
dt = scm.dt;

% initialize random number generator (from input argument)
rand_state = abs(1e3*scm.sim_num + round(scm.t0));
rng(rand_state, 'v5normal');

% number of time-steps for simulation
Nsteps = round(time_end / dt);
time = ( 0 : Nsteps - 1 )' * dt + scm.t0;

% noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
B_ee = scm.noise_sf .* sqrt(scm.noise_sc .* scm.phi_ee_sc / dt);
B_ei = scm.noise_sf .* sqrt(scm.noise_sc .* scm.phi_ei_sc / dt);

% % wave equation coefficients (as implemented in original Steyn-Ross sim)
% VLdt = SS.v * SS.Lambda * dt;
% c1 = VLdt.^2;
% c2 = 2 - 2*VLdt - VLdt.^2;
% c3 = 2*VLdt - 1;
% c4 = (SS.v * dt / dx).^2;

%% Define output variables

% Some are visualized, some are returned.

out_vars = scm.out_vars;

% Indices to capture from larger grid for NP and EC
indsNP = get_inds_(scm.grid_size, scm.centerNP, scm.dimsNP, 1);
indsEC = get_inds_(scm.grid_size, scm.centerEC, scm.dimsEC, scm.scaleEC);

% Initialize electrode output variables (structs NP, EC)
get_electrode_values;

%% Define dynamic variables
%Use as initial conditions the "last" values of previous simulation.

last = IC;  % initialize
new = last;  % initialize


%% Simulation
% new = last;
for ii = 1:Nsteps
    
	
    % Set the "source" locations' excitatory population resting voltage and
	% expand the wavefront
    update_drive
    
	% Update equations (update <new> values)
	update_wave_equations;
	update_synaptic_flux_eq;
	update_soma_voltages;
	update_firing_rates;
	update_extracellular_ion;
	update_gap_resting;
	
	% sanity check!
	if any(any(isnan(new.Vi)))
		error('Sigmoid generated NaNs!! (Either increase dx or reduce dt)');
	end
	
	% UPDATE the dynamic variables (pass <new> values to <last>).
	last = new;
    
    % Correct out-of-bounds values
	correct_OOB_values;
    
	get_electrode_values;
	visualize;
    
end

% Return requested variables.
no_return = out_vars(~ismember(out_vars, scm.return_fields));
NP = rmfield(NP, no_return);
EC = rmfield(EC, no_return);

% ------------------------------------------------------------------------
%% Nested dynamics functions

% 1. update wave equations
	function update_wave_equations
       
        new.phi2_ee = ...
			last.phi2_ee ...
			+ dt * ( ...
				-2 * scm.v .* scm.Lambda .* last.phi2_ee ...
				- (scm.v .* scm.Lambda).^2 .* last.phi_ee ...
				+ (scm.v .* scm.Lambda).^2 .* last.Qe...
			) ...
			+ dt * (scm.v / dx)^2 * del2_(last.phi_ee);
		 
		new.phi2_ei = last.phi2_ei ...
			+ dt * ( ...
				- 2*scm.v .* scm.Lambda .* last.phi2_ei ...
				- (scm.v .* scm.Lambda).^2 .* last.phi_ei ...
				+ (scm.v .* scm.Lambda).^2 .* last.Qe...
			) ...
			+ dt * (scm.v / dx)^2 * del2_(last.phi_ei);
							 
		new.phi_ee = last.phi_ee + dt * new.phi2_ee;
		new.phi_ei = last.phi_ei + dt * new.phi2_ei;
        
%         new.phi_ee = ...
%             c1*last.Qe + c2*last.phi_ee + c3*last.phi2_ee + ...
%             c4 * del2_(last.phi_ee);
%         new.phi2_ee = last.phi_ee;
%         
%         new.phi_ei = ...
%             c1*last.Qe + c2*last.phi_ei + c3*last.phi2_ei + ...
%             c4 * del2_(last.phi_ei);
%         new.phi2_ei = last.phi_ei;
		
	end

% 2. update the 4 synaptic flux equations (include sc noise)
	function update_synaptic_flux_eq

	%%%% E-to-E %%%%
		new.F_ee = last.F_ee ...
			+ dt * scm.gamma_e.^2 * ( ...
				- 2 / scm.gamma_e * last.F_ee ...
				- last.Phi_ee ...
				+ scm.Nee_a * last.phi_ee ...  % long range
				+ scm.Nee_b * last.Qe ...      % short range
				+ scm.noise_sc .* scm.phi_ee_sc ... % subcortical (tonic)
				+ B_ee .* randn(Nx, Ny) ...       % subcortical (random)
			);
		
		new.Phi_ee = last.Phi_ee + dt * new.F_ee;

	%%%% E-to-I %%%%
		new.F_ei = last.F_ei ...
			+ dt * scm.gamma_e.^2 * ( ...
				- 2 / scm.gamma_e * last.F_ei ...
				- last.Phi_ei ...
				+ scm.Nei_a * last.phi_ei ...    %long range
				+ scm.Nei_b * last.Qe ...   %short range
				+ scm.noise_sc .* scm.phi_ei_sc ...    %subcortical (tonic)
				+ B_ei .* randn(Nx, Ny)...   %subcortical (random)
			);
		
		new.Phi_ei = last.Phi_ei + dt * new.F_ei;

	%%%% I-to-E %%%%
		new.F_ie = last.F_ie ...
			+ dt * scm.gamma_i.^2 * ( ...
				- 2 / scm.gamma_i * last.F_ie ...
				- last.Phi_ie ...
				+ scm.Nie_b .* last.Qi ...     %short range
			);
		
		new.Phi_ie = last.Phi_ie + dt * new.F_ie;

	%%%% I-to-I %%%%
		new.F_ii = last.F_ii ...
			+ dt * scm.gamma_i.^2 * ( ...
				- 2 / scm.gamma_i * last.F_ii ...
				- last.Phi_ii ...
				+ scm.Nii_b .* last.Qi ...  %short range
			);
		new.Phi_ii = last.Phi_ii + dt * new.F_ii;

	end

% 3. update the soma voltages
	function update_soma_voltages
        
        if strcmpi(scm.drive_style, 'inhibitory')
            FS = scm.fixed_source(time(ii));
            IW = -scm.I_drive .* double(new.map * (time(ii) >=0));
        else  % excitatory drive style implements these in the dVe parameter
            [FS, IW] = deal(0);
        end

		new.Ve = last.Ve ...
			+ dt / scm.tau_e * ( ...
				(scm.Ve_rest - last.Ve) ...
				+ last.dVe + FS + IW ...                              % offsets from K+, FS, and IW
                + scm.rho_e * Psi_ee(last.Ve) .* last.Phi_ee ...      %E-to-E
                + scm.rho_i * Psi_ie(last.Ve) .* last.Phi_ie ...      %I-to-E
				+ last.Dee/dx^2 .* del2_(last.Ve) ...
			);
		new.Vi = last.Vi ...
			+ dt / scm.tau_i * ( ...
				(scm.Vi_rest - last.Vi) ...
				+ last.dVi + ...                                      % offset from K+
				+ scm.rho_e * Psi_ei(last.Vi) .* last.Phi_ei ...      %E-to-I
				+ scm.rho_i * Psi_ii(last.Vi) .* last.Phi_ii ...      %I-to-I
				+ last.Dii/dx^2 .* del2_(last.Vi) ...
			);
	end

% 4. update the firing rates
	function update_firing_rates
		
		new.Qe = ...
			scm.Qe_max .* ( ...
				1 ./ ( 1 + ...
					exp(-pi / ( sqrt(3) * scm.sigma_e ) .* ( last.Ve - scm.theta_e )) ...
				) ...
			) ; % No depo block on excitatory cells 
			
		new.Qi = ...
            scm.Qi_max_fun(new.state) .* ( ...
                1 ./ (1 + ...
                    exp(-pi / ( sqrt(3) * scm.sigma_i ) .* ( last.Vi - scm.theta_i )) ...
                ) ...
            ) ...     %The I voltage must be big enough,
            - scm.depo_block * scm.Qi_max_fun(new.state) .* ( ...
                1 ./ ( 1 + ...
                    exp(-pi / ( sqrt(3) * 3 ) .* ( last.Vi - ( scm.theta_i + 10 ) )) ...
                ) ...
            );   %... but not too big.

		
	end

% 5. Update extracellular ion.
	function update_extracellular_ion
        FR = last.Qe + last.Qi;
		new.K = ...
			last.K ...
			+ dt / scm.tau_K .* ( ...
				- scm.k_decay .* last.K ...  % decay term.
				+ scm.kR .* ...              % reaction term.
					FR ./ ( 1 + exp( -FR + 15) ) ...
                + scm.kD ./ dx^2 * del2_(last.K) ...  % diffusion term.
			);
        
        
        % This is not interacting with anything right now. Just thinking
        % about adding another ingredient to the sauce.
        new.GABA = ...  % chloride? 
			last.GABA ...
			+ dt / .1 .* ( ...
				- 1 .* last.GABA ...  % decay term.
				+ 1.5 .* ...              % reaction term.
					last.Qi ./ ( 1 + exp( -last.Qi ) ) ...
                + 1.5 * new.map ...
                + .25 ./ dx^2 * del2_(last.GABA) ...  % diffusion term.
			);
        
%         new.GABA = new.GABA + dt/2 * ( ...
%             1.5 * new.map ...  
%             - 1 * new.GABA ...
%             + .25/dx^2 * del2_(new.GABA));
%     
%         new.GABA = min(new.GABA, 1);
        
%         if all(size(SS.Nie_b) == params.grid_size)
%             new.GABA = SS.Nie_b;
%         end
			
	end

% 6. Update inhibitory gap junction strength, and resting voltages.  
   
	function update_gap_resting
       
        switch scm.gap_resting_update
            case 'dynamic'
                new.Dii = last.Dii + dt / scm.tau_dD * ( scm.KtoD .* last.K );
                new.Dee = new.Dii/100; 
                new.dVe = last.dVe + dt / scm.tau_dVe .* ( scm.KtoVe .* last.K );
                new.dVi = last.dVi + dt / scm.tau_dVi .* ( scm.KtoVi .* last.K );
            
            case 'sigmoid' % Dii ~ sigmoid(K); contant voltage offsets
                new.Dii = scm.sigmoid(last.K, scm.sigmoid_kD);
                new.Dee = new.Dii / 100;
                new.dVe = scm.sigmoid(last.K, scm.sigmoid_kdVe);
                new.dVi = scm.sigmoid(last.K, scm.sigmoid_kdVi);
            
        end
	end

% Correct out of bounds values
	function correct_OOB_values
		
		for f = fieldnames(last)'
			f = f{:}; %#ok<FXSET>
			if ismember(f, {'map', 'state'}) || all(isinf(scm.(f))), continue, end
			last.(f) = min(max(last.(f), scm.(f)(1)), scm.(f)(2));
		end
		
	end

% Update source and expand wavefront
	function update_drive
        
        switch scm.drive_style
            case {'', 'excitatory'}
                % Make sure I_drive and Qi_collapse are set to 0, nan or
                % else those will still show up in the voltage updates. 
                % Works with the Martinet FS parameters.
                if time(ii) > 0 
                    [new.map, new.state] = update_map_smooth( ...
                        last.state, scm.expansion_rate * dt / dx, ...
                        scm.excitability_map, dt);
                    if time(ii) <= scm.duration
                        new.dVe(new.map) = scm.source_drive; 
                        if isempty(scm.source), return, end
                        which_source = mod(floor(time(ii) / 2), size(scm.source, 3)) + 1;
                        new.dVe(scm.source(:, :, which_source)) = scm.source_drive;
                    else
                        if isnan(scm.post_ictal_source_drive), return; end			
                        new.dVe(new.map) = scm.post_ictal_source_drive;
                    end
                end
                
            case 'inhibitory'
                % The map for the IW updates here, but is applied directly
                % to the voltage dynamics rather than by changing dVi. This
                % makes it easier to modulate dVi with K+ if you want to do
                % that. Inhibitory collapse is applied here, however.
                if time(ii) > 0 
                    if isempty(scm.map)  % if no pregenerated map, do an update
                        [new.map, new.state] = update_map_smooth( ...
                            last.state, scm.expansion_rate * dt / dx, ...
                            scm.excitability_map, dt);
                    else
                        new.map = ...
                            scm.map > time(ii) - scm.excitability_map ...
                            & scm.map <= time(ii);
                        new.state = time(ii) - scm.map;
                    end
                    
                end
            otherwise
                error('Drive style ''%s'' not recognized');
        end
	end

%% Nested logistical functions

% Store electrode values
	function get_electrode_values
		
		if ~exist('NP', 'var')  % Initialize
			for v = union(out_vars, {'Qe', 'Ve'})
				NP.(v{:}) = zeros([Nsteps, scm.dimsNP], 'single');  % Microscale
				EC.(v{:}) = zeros([Nsteps, scm.dimsEC], 'single');  % Macroscale
			end
		else  % ... or store
			for v = union(out_vars, {'Qe', 'Ve'})
				NP.(v{:})(ii, :, :) = reshape(last.(v{:})(indsNP), [1, scm.dimsNP]); 
				
				% Take the local mean for ECOG electrodes
				blurEC = conv2(last.(v{:}), ones(3) ./ 9, 'same');  
				EC.(v{:})(ii, :, :) = reshape(blurEC(indsEC), [1, scm.dimsEC]);
			end
		end
	end

% Visualize results
	function visualize
		if scm.visualization_rate > 0
			if ~exist('fig', 'var') || isempty(fig)
				fig = create_fig_(out_vars, scm.grid_size, indsNP, indsEC);
                ss = mkdir(sprintf('SCM/%s/vids/%s_%d/', scm.label, scm.label, scm.sim_num));
			end
			if diff(floor((time(ii) - [dt 0]) * scm.visualization_rate))
				for sp = 1:numel(fig.ih)
					set(fig.ih(sp), 'cdata', last.(out_vars{sp}))
					colorbar
				end

				set(fig.ah, 'string', sprintf('T = %0.3f', time(ii)));
				drawnow;
                im = frame2im(getframe(fig.fig));
                imwrite(im, sprintf('SCM/%s/vids/%s_%d/%0.4f.png', scm.label, scm.label, scm.sim_num, time(ii)));
                
			end
		end
	end

%% Nested weighting functions

% e-to-e reversal-potential weighting function
	function weight = Psi_ee(V)
		weight = (scm.Ve_rev - V)./abs(scm.Ve_rev - scm.Ve_rest);
	end

% e-to-i reversal-potential weighting function
	function weight = Psi_ei(V)
		weight = (scm.Ve_rev - V)./abs(scm.Ve_rev - scm.Vi_rest);
	end

% i-to-e reversal-potential weighting function
	function weight = Psi_ie(V)
		weight = -(scm.Vi_rev - V)./max(abs(scm.Vi_rev - scm.Ve_rest), 1);
	end

% i-to-i reversal potential weighting function
	function weight = Psi_ii(V)
		weight = -(scm.Vi_rev - V)./max(abs(scm.Vi_rev - scm.Vi_rest), 1);
	end

end

% ------------------------------------------------------------------------
%% Supplementary functions


function [offset_i, offset_e] = dVi(t, gs, dVe)
    % Hack-y circular ring type IW
    if t < 0, offset_i = zeros(gs); offset_e = dVe; return; end
    WIDTH = 1;  % I think WIDTH probably needs to be at least Lambda to prevent excitatory transmission from inner to outer
    RATE = 1;  % in dx/s (so if dx = .1 cm, this is in mm/s)
    iw0 = [20 30];
    
    [xx, yy] = ndgrid(1:gs(1), 1:gs(2));
    xx = xx - iw0(1); yy = yy - iw0(2);
    D = sqrt(xx.^2 + yy.^2);
    
    r1 = D <= RATE * t & D >= RATE * t - WIDTH;  % outer ring
    r2 = D <= RATE * t - WIDTH & D >= RATE * t + (t/60 - 2) * WIDTH;  % inner ring  (40)
    r3 = D < RATE*t - 2*WIDTH;  % inside IW
    
    offset_i = zeros(gs);
    offset_i(r1) = 1000; % 5
    offset_i(r2) = 0; % -5
    offset_i(r3) = 0;
    offset_e = dVe;
%     offset_e(r1 | r2) = 2;
%     offset_e(r3) = 2;
    
    
    
end

%------------------------------------------------------------------------
function Y = del2_(X)

L = [0 1 0; 1 -4 1; 0 1 0];  % 5-point stencil Laplacian
% Y = convolve2(X, L, 'wrap');
% return
% r = 4/3 * (1 + 1 / sqrt(2));
% b = 1 / r; a = 1 / (sqrt(2) * r);
% a = .25; b = .5;
% L = [a b a; b -3 b; a b a];  % 9-point stencil Laplacian

% zero-flux BCs
	X = [repmat(X(1, :), 2, 1); X; repmat(X(end, :), 2, 1)];
	X = [repmat(X(:, 1), 1, 2), X, repmat(X(:, end), 1, 2)];

	Y = conv2(X, L, 'valid');
	Y = Y(2:end-1, 2:end-1);
end

%------------------------------------------------------------------------
function fig = create_fig_(out_vars, grid_size, addyNP, addyEC)
	stop_at = 'seizing_cortical_field>update_gap_resting';
	titles = out_vars;
	N = numel(titles);
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
		colormap(h(ii), 1-gray); axis equal; axis tight; axis ij; colorbar;
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
	legend(h(1), {'EC', 'NP'}, 'location', 'southeast', 'textcolor', [1 1 1])
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
function inds = get_inds_(grid_size, center, dims, scale)

	x_offset = ( (1:dims(1)) - floor(dims(1)/2) ) * scale;
	y_offset = ( (1:dims(2)) - floor(dims(2)/2) ) * scale;
	[xx, yy] = ndgrid(x_offset, y_offset);

	inds = sub2ind(grid_size, ...
		center(1) + xx(:), ...
		center(2) + yy(:));
end

%------------------------------------------------------------------------

