classdef SCM < handle
	
	methods  % constructor and other operators
		function scm = SCM(varargin)
			if nargin == 0, scm = scm.init(); return; end
			switch class(varargin{1})
				case 'SCM'
					scm = varargin{1};
				case 'struct'
					
					if isfield(varargin{1}, 'model')  % old params struct
						params = varargin{1};
					%	P.extract(params.meta);
						scm.IC = SCMState(params.IC);
						model = params.model;
						
						for ff = string(fieldnames(model)')
							if ~ismember(ff, properties(scm)), continue, end
							if isempty(model.(ff)), continue; end
							switch ff
								case "IC"
									continue
                                case {"noise" "sigmoids" "electrodes" "bounds" "time_constants"} 
									scm.extract(model.(ff));
								case "K"
									scm.potassium = model.K;
								otherwise
									scm.(ff) = model.(ff);
							end
							
						end
						scm = scm.extract(params.meta);
						
					else
						for arg = varargin
							scm = scm.extract(arg{:});
						end
					end
				case 'char'
                    if nargin == 1
                        mdl = validatestring(varargin{1}, ...
                            {'steyn-ross', 'martinet', 'wip', ...
                            'draft_Idriven', 'FS', 'SW', 'IW'});
                        switch mdl
                            case 'steyn-ross'
                                scm = scm.init();
                                scm.tau_dVe = Inf;  % voltage offsets are fixed
                                scm.tau_dVi = Inf;
                                scm.tau_dD = Inf;  % gap junction coupling is fixed
                                scm.depo_block = false;
                                scm.v = 140;
                                scm.IC.dVe = 1.5;
                                scm.IC.dVi = 0;
                                scm.D = .35;  % set this to whatever you like to change equilibrium
                                  % Note that lower than this doesn't
                                  % really look like the TW from the paper
                                  % though. Not sure what's wrong; could
                                  % just be that it needs firing injection
                                  % to push it into "up-state" (EDS.
                                  % 1/25/21)
                                scm.IC.Dii = scm.D;
                            case 'martinet'
                                % Fixed source Martinet model
                                scm = scm.init();
                                scm.sim_num = 1;
                                scm.padding = [5 10];
                                scm.out_vars = {'Qe', 'Ve', 'Qi', 'Vi', 'Dii', 'K'};
                                
                                scm.visualization_rate = 0;
                                scm.save = true;
                                scm.dx = .3;
                                scm.source_drive = 3;
                                
                                [xx, yy] = ndgrid(1:scm.grid_size(1), 1:scm.grid_size(2));
                                source = false(scm.grid_size);
                                source(abs(xx - scm.stim_center(1)) <= 1 & ...
                                    abs(yy - scm.stim_center(2)) <=1) = true;
                                scm.source = source;
                                
                                scm.expansion_rate = 0;  % No IW
                                scm.IC.dVe = 1;
                                
                                scm.noise_sf = 2;
                                scm.noise_sc = .2;
                            case 'FS'
                                % Generates a sim with a fixed source. Uses
                                % [10, 10] padding and a duration of 80
                                
                                scm = SCM('steyn');  % no dynamics on external drives
                                scm.label = 'FS';
                                scm.basename = 'SCM/FS/FS';
                                scm.dx = 0.1;
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);

                                scm.IC.dVi = 0;
                                scm.IC.dVe = 0;

                                dVe = 1;
                                scm.D = .3;  
                            
                                scm.save = true;
                                scm.visualization_rate = 0;
                                scm.depo_block = false;
                                scm.padding = [10 10];
                                scm.duration = 80;
                                
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.out_vars = {'Qe', 'Ve', 'dVe', 'K', 'Qi', 'Vi', 'dVi', 'Dii', 'map', 'state', 'GABA'};

                                % Design the IW
                                scm.expansion_rate = 0.1;  % 0.25
                                scm.excitability_map(scm.excitability_map > 0) = 3;
                                scm.I_drive = .3;


                                % Add a fixed source
                                center = round( [1.5 1.5] / scm.dx );
                                source_dims = round( [.2 .2] / scm.dx );

                                [xx, yy] = ndgrid(1:scm.grid_size(1), 1:scm.grid_size(2));
                                source = false(scm.grid_size);
                                source(abs(xx - center(1)) <= source_dims(1) & ...
                                    abs(yy - center(2)) <= source_dims(2)) = true;
                                scm.source = zeros(scm.grid_size);
                                scm.source(scm.excitability_map > 0) = dVe;
                                scm.source(source) = dVe + 1;


                                % Keep the electrodes in the center so that
                                % it's easy to rotate the source when
                                % generating sims for method comparison
                                scm.centerNP = round( scm.grid_size ./ 2 );


                                % Use Martinet dynamic potassium, but
                                % {dV*,Dii} follow sigmoid response
                                
            
                                % Add an IW source
                                scm.stim_center = round( [1.5 3.5] / scm.dx );  % 4.6 is edge

                                % Bound dVe/dVi
                                scm.dVe = [-Inf, 1];
                                scm.dVi = [-Inf, .3];

                                % No post ictal drive
                                scm.post_ictal_source_drive = nan;


                                scm.IC.Dii = scm.D;
                                scm.drive_style = 'inhibitory';
                                scm.map = scm.generate_map;
                                
                                
                                % Only simulate FS portion (i.e. update the
                                % IC to look like after IW passage)
                                % This was pretty interesting when running
                                % with dV*,Dii dynamic functions. Didn't do
                                % what I expected. (early waves had limited
                                % distance; sim cut off after 30(?)
                                % seconds; analyses did not detect TW after
                                % about 10 seconds)
%                                 scm.IC.Dii = .2;
%                                 scm.IC.dVe = .6;
%                                 scm.IC.dVi = .3;
%                                 scm.IC.K = .6;
%                                 scm.expansion_rate = 0;

                            case 'SW'
                                scm = SCM('FS');
                                scm.label = 'SW';
                                scm.basename = 'SCM/SW/SW';
                                temp = scm.source;
                                scm.rotate(pi/2);
                                scm.source = cat(3, scm.source, temp);


                            case 'IW'
                                % This is the sim that shows up in the last
                                % figure
                                
                                scm = SCM('steyn');  % no dynamics on external drives
                                scm.label = 'IW';
                                scm.basename = 'SCM/IW/IW';
                                scm.dx = 0.1;
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);
                                scm.sim_num = 0;
                                scm.IC.dVi = 0;
                                scm.IC.dVe = 0;
                                dVe = 1;
                                scm.save = true;
                                scm.visualization_rate = 0;
                                scm.depo_block = true;
                                scm.padding = [10 10];
                                scm.duration = 60;
                                
                                
                                % scm.dt = 1e-4;
                                
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.out_vars = {'Qe', 'Ve', 'dVe', 'K', 'Qi', 'Vi', 'dVi', 'Dii', 'map', 'state', 'GABA'};

                                % Design the IW
                                scm.expansion_rate = 0.1;  % 0.25
                                scm.excitability_map(scm.excitability_map > 0) = 3;
                                scm.I_drive = .3;
                                % scm.Nii_b = 100;


                                % Add a fixed source
                                center = round( [2.2 2.2] / scm.dx );
                                source_dims = round( [.2 .2] / scm.dx );

                                [xx, yy] = ndgrid(1:scm.grid_size(1), 1:scm.grid_size(2));
                                source = false(scm.grid_size);
                                source(abs(xx - center(1)) <= source_dims(1) & ...
                                    abs(yy - center(2)) <= source_dims(2)) = true;
                                scm.source = zeros(scm.grid_size);
                                scm.source(scm.excitability_map > 0) = dVe;
                                scm.source(source) = dVe + 1;


                                % Move the electrodes off the center
                                scm.centerNP = round( [3.5 3.5] / scm.dx );


                                % Potassium is dynamic but dV*, Dii have sigmoid response functions            

                                scm.stim_center = round( [2.0 3.5] / scm.dx );  % 4.6 is edge

                                scm.dVe = [-Inf, Inf]; % limits are naturally imposed with sigmoid functions
                                scm.dVi = [-Inf, Inf];

                                scm.post_ictal_source_drive = nan;

                                scm.drive_style = 'inhibitory';
                                scm.map = scm.generate_map;


                            case 'wip'
                                % Looking for a dVe/dVi pair that
                                % generates diffuse TW
                                
                                scm = SCM('steyn');  % no dynamics on external drives
                                scm.dx = 0.1;
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);
scm.sim_num = 9;
scm.t0_start = -2;
scm.IC.dVi = 0;
scm.IC.dVe = 0;
% scm.IC.dVe = randn(scm.grid_size) * .05;
% scm.IC.dVe = conv2(scm.IC.dVe, gausswin(10) * gausswin(10)', 'same');
dVe = 1;
scm.D = .35;  
                            scm.save = true;
                            scm.visualization_rate = 10;
                                scm.depo_block = true;
                                scm.padding = [10 10];
scm.duration = 60;
                                
                                
scm.dimsNP = [4 4];
                                % scm.dt = 1e-4;
                                
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.out_vars = {'Qe', 'Ve', 'dVe', 'K', 'Qi', 'Vi', 'dVi', 'Dii', 'map', 'state', 'GABA'};

                                % Design the IW
                                scm.expansion_rate = 0.1;  % 0.25
% scm.excitability_map(scm.excitability_map > 0) = 3;
scm.I_drive = .3;
% scm.Nii_b = 100;


% Add a fixed source
center = round( [2.2 2.2] / scm.dx );
source_dims = round( [.2 .2] / scm.dx );

[xx, yy] = ndgrid(1:scm.grid_size(1), 1:scm.grid_size(2));
source = false(scm.grid_size);
source(abs(xx - center(1)) <= source_dims(1) & ...
    abs(yy - center(2)) <= source_dims(2)) = true;
scm.source = zeros(scm.grid_size);
scm.source(scm.excitability_map > 0) = dVe;
scm.source(source) = dVe + 1;


% Move the electrodes off the center
scm.centerNP = round( [3.5 3.5] / scm.dx );


% Add Martinet Potassium dynamics (commented out because testing sigmoid
% responses instead)
% % scm.IC.K = .3;
% scm.tau_dVe = 250 * 1;  % Speed up Ve reaction (with 0.8)
% scm.tau_dVi = 250 * 1;  % Speed up Vi reaction (with 0.8)
% scm.tau_dD = 200 * 50;  % slow down the changes

                              

scm.stim_center = round( [2.0 3.5] / scm.dx );  % 4.6 is edge

% scm.dVe = [-Inf, .7];  % built into sigmoids
% scm.dVi = [-Inf, .3];



% scm.IC.dVe = zeros(scm.grid_size);
% scm.IC.dVe(scm.excitability_map > 0) = dVe;
scm.post_ictal_source_drive = nan;



% scm.Dii = [0.2 Inf];
scm.IC.Dii = scm.D;  % not currently dynamic... function of K
scm.drive_style = 'inhibitory';
scm.excitability_map = 3 * ones(scm.grid_size);
% scm.IC.Dii = ndgrid(linspace(.01, .4, 50), 1:50);
% scm.IC.Dee = scm.D/100;  % not used in Martinet formulation (always
% updates to Dii/100)

% scm.IC.Qe = 18.47;
% scm.IC.Qi = 32.68;
% scm.IC.Ve = -57.71;
% scm.IC.Vi = -58.01;


% scm.map = scm.generate_map;

                                
                                
                            otherwise
                                error('Input ''%s'' not recognized.', mdl)
                        end
                    else
                        for ii = 1:2:nargin
                            ff = validatestring(varargin{ii}, fieldnames(scm));
                            scm.(ff) = varargin{ii + 1};
                        end
                    end
			end
            scm = scm.init();
			
		end
        function P = init(P)
            if isempty(P.basename)  % set default basename
				str = sprintf('%s/%s/%s', P.label, P.label, P.label); 
				P.basename = str;
            end
%             if ~isstruct(P.seed)  % set seed
% 				P.seed = rng(P.seed);
%             end
            if isempty(P.IC), P.IC = []; end  % set default IC
			P.IC = P.IC.resize(P.grid_size);
			if all(P.stim_center) && isempty(P.map)
				P.IC.map = false(P.grid_size);  % Source of ictal activity (on/off)
                if all(P.stim_center > 0)
                    P.IC.map(P.stim_center(1), P.stim_center(2)) = 1;
                end
                P.IC.state = double(P.IC.map) - 1;  % Seizure state (ictal/non-ictal)
            elseif ~isempty(P.map)
                P.IC.map = ...
                    P.map > P.t0_start - P.excitability_map ...
                    & P.map <= P.t0_start;
                P.IC.state = P.t0_start - P.map;
            end
            
        end

        
        function clean(scm)
            delete(sprintf('%s_%d_*.mat', scm.basename, scm.sim_num));
        end
        function map_ = generate_map(scm)
            % Generates a linear radial distance from source map with some
            % 2D gaussian noise
            [xx, yy] = find(ones(scm.grid_size));
            xx = xx - scm.stim_center(1);
            yy = yy - scm.stim_center(2);

            T = reshape( ...
                scm.expansion_rate / scm.dx * sqrt(xx.^2 + yy.^2), ...
                scm.grid_size);

            noise_ = randn(size(T)) * 2;
            f_smooth = round(.5/scm.dx);
            win = gausswin(f_smooth) * gausswin(f_smooth)';
            map_ = T + conv2(noise_, win./sum(win, 'all'), 'same');
            map_(map_ < 0) = 0;
%             map_(~scm.excitability_map) = Inf;

        end
        function rotate(self, theta), self.Rotate(theta); end  % legacy
        function inds = NPinds(scm)

            x_offset = ( (1:scm.dimsNP(1)) - floor(scm.dimsNP(1)/2) );
            y_offset = ( (1:scm.dimsNP(2)) - floor(scm.dimsNP(2)/2) );
            [xx, yy] = ndgrid(x_offset, y_offset);

            inds = sub2ind(scm.grid_size, ...
                scm.centerNP(1) + xx(:), ...
                scm.centerNP(2) + yy(:));
        end
        function inds = ECinds(scm)
            % Never used this... might need help
            x_offset = ( (1:scm.dimsNP(1)) - floor(scm.dimsNP(1)/2) ) * scm.scaleEC;
            y_offset = ( (1:scm.dimsNP(2)) - floor(scm.dimsNP(2)/2) ) * scm.scaleEC;
            [xx, yy] = ndgrid(x_offset, y_offset);

            inds = sub2ind(scm.grid_size, ...
                scm.centerNP(1) + xx(:), ...
                scm.centerNP(2) + yy(:));
        end
        function show_layout(scm)
            X = zeros(scm.grid_size);
            
            X(scm.stim_center(1), scm.stim_center(2)) = 1;
            X(scm.NPinds) = 2;
            for ii = 1:size(scm.source, 3)
                temp = scm.source(:, :, ii);
                X(temp == max(temp, [], 'all')) = ii + 2;
            end
            X(~scm.excitability_map) = -1;
            figure; imagesc(X);
            cb = colorbar;
            set(cb, 'ticks', 0:3, 'limits', [-1 ii+2], ...
                'ticklabels', {'EZ', 'IWs', 'MEA', 'FS'})
        end
        function phi = theta_FS(scm)
            phi = nan(size(scm.source, 3), 1);
            for ii = 1:numel(phi)
                [xx, yy] = find(scm.source(:, :, ii) == max(scm.source(:, :, ii), [], 'all'));
                dy = scm.centerNP(2) - mean(yy);
                dx_ = scm.centerNP(1) - mean(xx);
                phi(ii) = atan2(dy, dx_);
            end
        end
        function phi = theta_IW(scm)
            dy = scm.centerNP(2) - scm.stim_center(2);
            dx_ = scm.centerNP(1) - scm.stim_center(2);
            phi = atan2(dy, dx_);
        end
        function Run(self)
            self.CreateDirectory;
            self.RunSimulation;
            self.ConvertToMea;
        end
        

		Rotate(self, theta)
        CreateDirectory(self)
        RunSimulation(self)
        ConvertToMea(self)
        AnalyzeWaveDirections(self)
        Preview(self)
	end
    
    methods (Access = private)
        map = DefaultExcitabilityMap(P)
    end
    properties (Access = private)
        Qe_movie
        Ve_movie
        K_movie
    end
	
	properties  
        
        
        map  % a predetermined map. Generate this if the IW is independent of the rest of the sim
        
        % meta
		basename
		sim_num
		save = true  % Save output
		visualization_rate = 0  % Show this many frames per second
		t_step = 1  % Simulate <t_step> second intervals
		t0_start  % Allows continue from previous sim (IC will use last from file number t0_start-1)
		duration (1,1) double {mustBeNonnegative} = 60
		padding (1,2) = [10 10]  % Padding before and after seizure
		source_drive (1, 1) double = 2.5
		post_ictal_source_drive (1,1) double = 1.5
        I_drive (1, 1) double = .7  % Allow inhibitory driving (EDS, 1/26/21)
        post_ictal_I_drive (1, 1) double = nan  % Allow inhibitory driving (EDS, 1/26/21)
		subsample (1,1) double {mustBePositive} = Inf  % Allow downsampling when creating mea
		return_fields (1,:) = {'Qe', 'Ve'}  % Qe required to make mea
		out_vars (1,:) = {'Qe', 'Ve'}  % Define which variables you would like to visualize (can be any in IC)
		source % Define fixed, alternating sources of excitation
% 		seed = rng  % Set seed for repeatability
		label (1,:) char = 'SCM'
		t0  % used to keep track of progress in sims 
        drive_style (1,:) char = 'excitatory'  % Allow inhibitory driving (EDS, 1/26/21)
        
        % model
		stim_center (1,2) = [20 20]
		grid_size (1,2) {mustBePositive} = [50 50] % size of grid to simulate (must be even)
		dt = 2e-4
		dx = .4  % (cm) 
		expansion_rate (1,1) double {mustBeNonnegative} = .625  % in cm^2/s; set to 0 for fixed source
		excitability_map  % boundary to IW spread
		IC SCMState 
        
        
        % steady states
		tau_e = 0.04  % excit neuron time-constant (/s) [original = 0.04 s, Martinet = 0.02 s]
		tau_i = 0.04  % inhib neuron time-constant (/s) [original = 0.04 s, Martinet = 0.02 s]

		% voltage limits
		Ve_rev = 0  % reversal potential (mV)
		Vi_rev = -70
		Ve_rest = -64  % resting potential (mV)
		Vi_rest = -64

		% gain per synapse at resting voltage (millivolt.sec)
		rho_e = 1.00e-3
		rho_i = -1.05e-3

		% E/IPSP rate constants
		gamma_e = 170  % EPSP decay rate (/s)
		gamma_i = 50  % IPSP decay rate (/s)

		% gap-junction diffusive-coupling strength
		D = 0.8  % i <--> i steady state (cm^2) [0, 1]  [Martinet = 1]

		% sigmoid characteristics
		Qe_max = 30  % sigmoid maximum (s^-1)
		Qi_max = 60
		theta_e = -58.5	 % sigmoid threshold (mV)
		theta_i = -58.5
		sigma_e = 3.0	% sigmoid 'width' (mV)
		sigma_i = 5.0
        depo_block = 1  % apply depolarization block

		% connectivities: j-->k convention (dimensionless)			
		Nee_a = 2000  % cortico-cortical
		Nei_a = 2000
		Nee_b = 800
		Nei_b = 800
		Nie_b = 600
		Nii_b = 600
        Nee_sc = 50
        Nei_sc = 50
		

		% axonal conduction velocity (cm/s), 
		v = 280  % [original = 140 cm/s]

		% inverse-length scale for connectivity (/cm)
		Lambda = 4.0	
        
        
        % time_constants
		tau_dD = 200  %inhibitory gap junction time-constant (/s).
		tau_dVe = 250  %excitatory population resting voltage time-constant (/s).
		tau_dVi = 250  %inhibitory population resting voltage time-constant (/s).
        
        
        % potassium
		tau_K = 200    %time-constant (/s).
		k_decay = 0.1  %decay rate (/s).
		kD = 1       %diffusion coefficient (cm^2/s).
		KtoVe = 10     %impact on excitatory population resting voltage.
		KtoVi = 10     %impact on inhibitory population resting voltage.
		KtoD = -50    %impact on inhibitory gap junction strength.
		kR = 0.15   %scale reaction term. 

		% Potassium excitability function
		k_peak = .5  % K concentration at which voltage offset is greatest
		k_width = .1  % Width of K concentration peak

		
        % electrodes
		centerNP  % defaults to grid center
		dimsNP = [10 10]

		% Macroscale
		centerEC  % defaults to grid center
		scaleEC = 4
		dimsEC = [3 3]
        
        
        % bounds
		Dee (1,2) double = [-Inf Inf]  % i <--> i gap-junction diffusive-coupling strength (cm^2)
		Dii (1,2) double = [.009 Inf]  % The inhibitory gap junctions cannot pass below a minimum value of 0.009 / dx^2.
		K (1,2) double = [-Inf 1] % extracellular potassium concentration (cm^2)
        GABA (1, 2) double = [-Inf Inf]  % GABA
		Qe (1,2) double = [0 Inf]  % Activity of excitatory population.
		Qi (1,2) double = [0 Inf]  % Activity of inhibitory population.
		Ve (1,2) double = [-Inf Inf]  % Voltage  of excitatory population.
		Vi (1,2) double = [-Inf Inf]  % Voltage of inhibitory population.
		dVe (1,2) double = [-Inf 1.5]  % Excitatory resting potential offset (mV)
		dVi (1,2) double = [-Inf 0.8]  % Inhibitory resting potential offset (mV)
		Phi_ee (1,2) double = [-Inf Inf]  % e <--> e synaptic flux
		Phi_ei (1,2) double = [-Inf Inf]  % e <--> i synaptic flux
		Phi_ie (1,2) double = [-Inf Inf] % i <--> e synaptic flux
		Phi_ii (1,2) double = [-Inf Inf]  % i <--> i synaptic flux
		phi2_ee (1,2) double = [-Inf Inf]  % Wave dynamics
		phi2_ei (1,2) double = [-Inf Inf]
		phi_ee (1,2) double = [-Inf Inf]
		phi_ei (1,2) double = [-Inf Inf]
		F_ee (1,2) double = [-Inf Inf]  % flux dynamics
		F_ei (1,2) double = [-Inf Inf]
		F_ie (1,2) double = [-Inf Inf]
		F_ii (1,2) double = [-Inf Inf]

        
        
        % noise
		noise_sf (1,1) double {mustBeNonnegative} = 4  % [Martinet = 2]
		noise_sc (1,1) double {mustBeNonnegative} = 0.2;
        
        
        % sigmoids
        % Don't use this style anymore. Use sigmoid_xxx instead
        % i.e. sigmoid_xxx = [xxx_center params.xxx(2) xxx_width*4]
        
% 		kdVe_center = 0.8  % center of K-->dVe sigmoid  [Martinet = 0.8]
% 		kdVe_width = .1  % ... and width  [Martinet = 0.8]
% 		kD_center = 3e-4  % center of K-->Dii sigmoid [Martinet = .06 / .85]
% 		kD_width =  5e-4  % ... and width [Martinet = .01 / .06]
        
        sigmoid_kdVe = [0.5 .7 10]  % [mid max slope] (slope is ~max/width*4 if you prefer to think of it that way)
        sigmoid_kdVi = [.5 .3 10]
        sigmoid_kD = [1 .3 -4]

	end

	
	properties (Dependent = true)
        
        % default subcortical fluxes
        % %% ORIGINAL FORMULATION %%
        % [Nee_sc,Nei_sc]= deal(50, 50)  % subcortical  
		% phi_ee_sc = Nee_sc * Qe_max  % original [1500]
		% phi_ei_sc = Nei_sc * Qe_max  % original [1500]

        % %% EDS %%
		phi_ee_sc
		phi_ei_sc
        
		% d/dV derivatives of psi_ij weighting functions
		d_psi_ee 
		d_psi_ei 
		d_psi_ie 
		d_psi_ii 
		
		% Nee and Nie totals for cortico-cortical plus intracortical
		Nee_ab
		Nei_ab
        
        
        
        % These are here to match old parameter structures, but I think
        % I've stopped keeping up with legacy so probably pointless
        meta
		model
		electrodes
		time_constants
		potassium
		noise
		sigmoids
		bounds
		SS


    end
    
    
    
	methods  % model parameter getters
        function phi_ee_sc = get.phi_ee_sc(self)
            phi_ee_sc = self.Nee_sc * self.Qe_max;
        end
        function phi_ei_sc = get.phi_ei_sc(self)
            phi_ei_sc = self.Nei_sc * self.Qe_max;
        end

		function d_psi_ee = get.d_psi_ee(S)
			d_psi_ee = -1./(S.Ve_rev - S.Ve_rest);
		end
		function d_psi_ei = get.d_psi_ei(S)
			d_psi_ei = -1./(S.Ve_rev - S.Vi_rest);
		end
		function d_psi_ie = get.d_psi_ie(S)
			d_psi_ie = -1./(S.Vi_rev - S.Ve_rest);
		end
		function d_psi_ii = get.d_psi_ii(S)
			d_psi_ii = -1./(S.Vi_rev - S.Vi_rest);
		end
		
		% Nee and Nie totals for cortico-cortical plus intracortical
		function N = get.Nee_ab(S)
            N = S.Nee_a + S.Nee_b;
		end
		function N = get.Nei_ab(S)
            N = S.Nei_a + S.Nei_b;
		end

    end

    methods  % getters for meta
		
        function t0s = get.t0_start(p)
            if isempty(p.t0_start), p.t0_start = -p.padding(1); end
            t0s = p.t0_start;
        end
		
		
		
		function num = get.sim_num(P)
			s = 1;
			while isempty(P.sim_num)
				if isempty(dir([P.basename '_' num2str(s) '*.mat']))
					P.sim_num = s; 
				end
				s = s+1;
			end
			num = P.sim_num;
        end
        
        
        % model

		function dt = get.dt(P)
			if P.dt > 2e-4 && (any(P.IC.Dii(:) >= 0.87) || P.SS.v(:) > 140)
				P.dt = 2e-4;
				warning('High Dii or SS.v; setting dt to 2e-4.');
			end
			dt = P.dt;

		end
		function grid_size = get.grid_size(p)
			if any(mod(p.grid_size, 2)), p.grid_size = p.grid_size + mod(p.grid_size, 2); end
			grid_size = p.grid_size;
		end
		function map = get.excitability_map(P)
			if isempty(P.excitability_map)
				P.excitability_map = P.DefaultExcitabilityMap; 
			end
            map = P.excitability_map;

			assert(all(size(map) == P.grid_size))
		end
		
	end
	

	methods  % electrodes
		function center = get.centerNP(P)
            if isempty(P.centerNP)
				P.centerNP = round(P.grid_size / 2); 
            end
            center = P.centerNP;
		end
		function center = get.centerEC(P)
            if isempty(P.centerEC)
                P.centerEC = round(P.grid_size / 2);
            end
            center = P.centerEC;
		end
		
	end
	
	
	methods  % dependent getters
		function obj = substruct(P, names)
			for ff = names
				obj.(ff) = P.(ff);
			end
		end
		function P = extract(P, obj)
			for ff = string(fieldnames(obj)')
				P.(ff) = obj.(ff);
			end
		end
		function set.meta(P, value), P.extract(value); end
		function set.potassium(P, value), P.extract(value); end
		function set.noise(P, value), P.extract(value); end
		function set.model(P, value), P.extract(value); end
		function set.electrodes(P, value), P.extract(value); end
		function set.bounds(P, value), P.extract(value); end
		function set.sigmoids(P, value), P.extract(value); end
		function set.SS(P, value)
			for ff = string(fieldnames(value)')
				if strcmp(ff, 'Dii'), P.D = value.Dii;
				elseif strcmp(ff, 'Dee'), continue;
				else, P.(ff) = value.(ff);
				end
			end
		end
		function meta = get.meta(P)
			names = ["basename", "sim_num", "save", "visualization_rate", ...
                "t_step", "t0_start", "duration", "padding", ....
                "source_drive", "post_ictal_source_drive", "subsample", ...
                "return_fields", "out_vars", "source", "label", "drive_style"];  % "seed", 
			meta = P.substruct(names);
		end
		function model = get.model(P)
			names = ["stim_center", "grid_size", "dt", "dx", ...
                "expansion_rate", "excitability_map"];
			model = P.substruct(names);
		end
		function time_constants = get.time_constants(P)
			names = ["tau_dD", "tau_dVe", "tau_dVi"];
			time_constants = P.substruct(names);
		end
		function potassium = get.potassium(P)
			names = ["tau_K", "k_decay", "kD", "KtoVe", "KtoVi", "KtoD", ...
                "kR", "k_peak", "k_width"];
			potassium = P.substruct(names);
		end
		function electrodes = get.electrodes(P)
			names = ["centerNP", "dimsNP", "centerEC", "scaleEC", "dimsEC"];
			electrodes = P.substruct(names);
		end
		function bounds = get.bounds(P)
			names = ["Dee", "Dii", "K", "Qe", "Qi", "Ve", "Vi", "dVe", ...
                "dVi", "Phi_ee", "Phi_ei", "Phi_ie", "Phi_ii", ...
                "phi2_ee", "phi2_ei", "phi_ee", "phi_ei", "F_ee", ...
                "F_ei", "F_ie", "F_ii", "GABA"];
			bounds = P.substruct(names);
		end
		function noise = get.noise(P)
			names = ["noise_sf", "noise_sc"];
			noise = P.substruct(names);
		end
		
		function SS = get.SS(P)
			names = ["tau_e", "tau_i", "Ve_rev", "Vi_rev", "Ve_rest", ...
                "Vi_rest", "rho_e", "rho_i", "gamma_e", "gamma_i", "D", ...
                "Qe_max", "Qi_max", "theta_e", "theta_i", "sigma_e", ...
                "sigma_i", "Nee_a", "Nei_a", "Nee_b", "Nei_b", "Nie_b", ...
                "Nii_b", "phi_ee_sc", "phi_ei_sc", "v", "Lambda", ...
                "d_psi_ee", "d_psi_ei", "d_psi_ie", "d_psi_ii", ...
                "Nee_ab", "Nei_ab", "depo_block"];
			SS = P.substruct(names);
			SS.Dii = SS.D;  % naming used in seizing_cortical_field
		end
	end
	
	
	%% Static methods
	methods (Static)
		
		function center = get_center(center, dims, grid_size)
			% (Static) center = get_center(center, dims, grid_size) 
			if isempty(center)
				center = round(grid_size / 2); 
			end
			center = min(max(center, ceil(dims / 2)), grid_size - floor(dims / 2));
		end
	end
	

end

%%

% 		function ic = get.IC(P)
% 			%% Initial conditions
% 			% These are the same as in Waikato-Kramer except that dVe defaults to -1
% 			% instead of 1 (original IC in <map_type='ictal_wavefront'> scenario).
% 			% Additionally, resizing is done via bootstrapping in the case of a grid
% 			% size different from 100x100.
% 
% 			
% 			ic.Dee = P.SS.Dee;  % i <--> i gap-junction diffusive-coupling strength (electrodes)
% 			ic.Dii = P.SS.Dii;  % e <--> e gap-junction diffusive-coupling strength in all space (electrodes)
% 			ic.K = 0;  % extracellular potassium concentration (cm^2)
% 			ic.Qe = 0;  % Activity of excitatory population.
% 			ic.Qi = 0;  % Activity of inhibitory population.
% 			ic.Ve = P.SS.Ve_rest;  % Voltage  of excitatory population.
% 			ic.Vi = P.SS.Vi_rest;  % Voltage of inhibitory population.
% 			ic.dVe = -1;  % Excitatory resting potential offset (mV)
% 			ic.dVi = .1;  % Inhibitory resting potential offset (mV)
% 			ic.Phi_ee = 0;  % e <--> e synaptic flux
% 			ic.Phi_ei = 0;  % e <--> i synaptic flux
% 			ic.Phi_ie = 0;  % i <--> e synaptic flux
% 			ic.Phi_ii = 0;  % i <--> i synaptic flux
% 			ic.phi2_ee = 0;  % Wave dynamics
% 			ic.phi2_ei = 0;
% 			ic.phi_ee = 0;
% 			ic.phi_ei = 0;
% 			ic.F_ee = 0;  % flux dynamics
% 			ic.F_ei = 0;
% 			ic.F_ie = 0;
% 			ic.F_ii = 0;
% 			
% 			for ff = string(fieldnames(ic)')
% 				ic.(ff) = ic.(ff) .* ones(P.grid_size); 
% 			end
% 			
% 			ic.map = false(P.grid_size);  % Source of ictal activity (on/off)
% 			ic.map(P.stim_center(1), P.stim_center(2)) = 1;
% 			ic.state = double(ic.map);  % Seizure state (ictal/non-ictal)
% 			
% 
% 		end

