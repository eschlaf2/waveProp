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
                            'FS', 'SW', 'IW', 'IW_c7', 'IW_mg63', 'IW_coloc', ...
                            'IW_big'});
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
                                scm.drive_style = 'excitatory';
                                scm.I_drive = 0;
                                scm.Qi_collapse = nan;
                                scm.K = [-Inf 1];
                                
                                scm.visualization_rate = 0;
                                scm.save = true;
                                scm.dx = .3;
                                scm.source_drive = 3;
                                scm.gap_resting_update = 'dynamic';
                                
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
                                
                                scm = SCM('IW');  % no dynamics on external drives
                                scm.label = 'FS';
                                scm.basename = 'SCM/FS/FS';
                                scm.dx = 0.1;
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);

                                
                                scm.save = true;
                                scm.visualization_rate = 0;
                                scm.depo_block = false;
                                scm.padding = [10 10];
                                scm.duration = 32;
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.out_vars = {'Qe', 'Ve', 'dVe', 'K', 'Qi', 'Vi', 'dVi', 'Dii', 'map', 'state', 'GABA'};

                                % Design the IW
%                                 scm.expansion_rate = 0.1;  % 0.25
%                                 scm.excitability_map(scm.excitability_map > 0) = 3;
%                                 scm.I_drive = .3;
                                scm.expansion_rate = 0;


                                % Add a fixed source. Use the same value as
                                % in 'IW', but change the location
                                center = round( [1.5 1.5] / scm.dx );
                                radius = round( [.25 .25] / scm.dx );

                                dVe = max(scm.source, 'all');  
                                scm.source = double(scm.ellipse(center, radius)) * dVe;


                                % Keep the electrodes in the center so that
                                % it's easy to rotate the source when
                                % generating sims for method comparison
                                scm.dimsNP = [4 4];
                                scm.centerNP = round( scm.grid_size ./ 2 );


%                                 % Add an IW source
%                                 scm.stim_center = round( [1.5 3.5] / scm.dx );  % 4.6 is edge

                                % No post ictal drive
                                scm.post_ictal_source_drive = nan;

                                
                                scm.drive_style = 'inhibitory';
                                
                                
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
                                scm.dimsNP = [4 4];
                                scm.dt = 2e-4;
                                scm.grid_size = round( [5 5] / scm.dx);
                                scm.sim_num = 0;
                                scm.save = true;
                                                               
                                scm.visualization_rate = 10;
                                scm.I_drive = 0.5;
                                scm.depo_block = false;
                                scm.padding = [10 10];
                                scm.duration = 60;                                
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.out_vars = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};

                                % Design the IW
                                scm.expansion_rate = 0.1;  % 0.25
                                scm.stim_center = round( [2.0 3.5] / scm.dx );  % 4.6 is edge
                                scm.excitability_map = 3 * ones(scm.grid_size);


                                % Add a fixed source and IZ
                                scm.IZ = 1.2;
                                center = round( [2.2 2.2] / scm.dx );
                                radius = round( [.25 .25] / scm.dx );
                                scm.source = scm.ellipse(center, radius, 1);
                                

                                % Move the electrodes off the center
                                scm.centerNP = round( [3.5 3.5] / scm.dx );


                                % Potassium is dynamic but dV*, Dii have sigmoid response functions            

                                

                                scm.dVe = [-Inf, Inf]; % limits are naturally imposed with sigmoid functions
                                scm.dVi = [-Inf, Inf];

                                scm.post_ictal_source_drive = nan;

                                scm.drive_style = 'inhibitory';

                            case 'IW_big'
                                scm = SCM('IW');
                                scm.label = 'IW_big';
                                scm.basename = 'SCM/IW_big/IW_big';
                                scm.dimsNP = [24 24];
                                scm.centerNP = scm.stim_center;
                                
                                scm.visualization_rate = 0;
                                scm.duration = 35;
                                scm.padding = [2 2];
                                
                                % make a small version to compare (upper
                                % left corner of big MEA)
%                                 scm.dimsNP = [4 4];
%                                 scm.centerNP = [9 24];
%                                 scm.visualization_rate = 10;
                                
                                
                            case 'IW_coloc'
                                scm = SCM('IW');
                                scm.label = 'IW_coloc';
                                scm.basename = 'SCM/IW_coloc/IW_coloc';
                                
                                scm.source = scm.ellipse(scm.stim_center, 3, 1);
                                scm.visualization_rate = 0;

                            case 'IW_c7'
                                
                                scm = SCM('IW');  % no dynamics on external drives
                                scm.sim_num = 90;
                                scm.I_drive = 0;
                                scm.set_layout('d');
                                
                                scm.centerNP = [26.5 6.5];
                                scm.dimsNP = [4 4];
                                scm.map = scm.generate_contagion;
                                
                            case 'IW_mg63'
                                
                                scm = SCM('IW');  % no dynamics on external drives
                                scm.sim_num = 91;
                                scm.I_drive = 0;
                                scm.set_layout('c');
                                
                                scm.centerNP = [26.5 30.5];
                                scm.dimsNP = [4 4];
                                scm.map = scm.generate_contagion;
                                
                            case 'wip'
                                % Look at stronger FS
                                
                                scm = SCM('IW');  % no dynamics on external drives
                                
scm.sim_num = 9;

                            scm.visualization_rate = 10;
                                scm.padding = [2 10];
scm.duration = 60;
scm.dimsNP = [4 4];                                           
                                
                                
                                
                                % Save and visualize some extra fields for testing
                                scm.return_fields = {'Qe', 'Ve', 'Qi', 'Vi', 'K', 'Dii'};
                                scm.out_vars = {'Qe', 'Ve', 'dVe', 'K', 'Qi', 'Vi', 'dVi', 'Dii', 'map', 'state', 'GABA'};

scm.stim_center = [25 20];
scm.centerNP = [25 15];
% scm.IZ = scm.ellipse([], [], 1.2);
scm.source = scm.ellipse([25 10], 2, .8);
scm.sigmoid_kD(2) = .2;
% scm.expansion_rate = 0;
% scm.excitability_map = 5 * ones(scm.grid_size);
% scm.Qi_collapse([1 2]) = [30 -3];
% scm.I_drive = 2;
scm.t0_start = 0;


% Add Martinet Potassium dynamics (commented out because testing sigmoid
% responses instead)
% % scm.IC.K = .3;
% scm.tau_dVe = 250 * 1;  % Speed up Ve reaction (with 0.8)
% scm.tau_dVi = 250 * 1;  % Speed up Vi reaction (with 0.8)
% scm.tau_dD = 200 * 50;  % slow down the changes

% default behavior is to do this as part of the integration; here, this is to
% override the map from the last time step (so you don't have to rerun the
% burn-in when you change the map)
% scm.map = scm.generate_contagion;  
scm.map = scm.generate_map; % This works ok with scm.Qi_collapse low
% (e.g. 10)

                                
                                
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

        function set_layout(scm, idx)
            
            if nargin < 2, idx = 'a'; end
            iz = scm.IZ;
            [xx, yy] = find(iz > 0);
            C = mean([xx yy]);  % IZ center
            R = max(abs([xx yy] - C));  % IZ radius
            fsr = .25/scm.dx;  % fixed source radius
            
            r = [.3 .1] .* R;
            switch lower(idx)
                case 'a'
                    scm.stim_center = C + [1 -1] .* r;
                    scm.source = scm.ellipse(C + [1 1] .* r, fsr, 1);
                case 'b'
                    scm.stim_center = C + [1 -10] .* r;
                    scm.source = scm.ellipse(C + [1 -8] .* r, fsr, 1);
                case 'c'
                    scm.stim_center = C + [1 -2] .* r;
                    scm.source = scm.ellipse(C + [1 6] .* r, fsr, 1);
                case 'd'
                    scm.stim_center = C + [1 -10] .* r;
                    scm.source = scm.ellipse(C + [1 -2] .* r, fsr, 1);
                otherwise
                    error('Input ''%s'' not recognized', idx)
            end
        end
        function set_random_layout(scm)
            rng(scm.sim_num+500); % set the seed
            rand(30);  % blow out some numbers
            
            dVe_base = mode(scm.source, 'all');
            dVe_max = max(scm.source, [], 'all');

            scm.source = double(scm.ellipse()) * dVe_base;
            iz_inds = find(scm.ellipse());  % define the IZ
            N = numel(iz_inds);
            gs = scm.grid_size;
            [xx, yy] = meshgrid(1:gs(1), 1:gs(2));

            iw_loc = iz_inds(randi(N));
            fs_loc = iz_inds(randi(N));
            mea_loc = iz_inds(randi(N));

            scm.stim_center = [xx(iw_loc), yy(iw_loc)];  % get IW onset location
            scm.source(scm.ellipse([xx(fs_loc), yy(fs_loc)], 3)) = dVe_max;
            scm.centerNP = [xx(mea_loc) yy(mea_loc)];
            
            scm.init;
                                
        end
        function fs = fixed_source(scm, t)
            % Get the voltage offset for the fixed source at time t
            % If there are multiple fixed sources, rotate through them
            % every W seconds.
            W = 2;
            which_source = mod(floor(t / W), size(scm.source, 3)) + 1;
            if t >=0 && t < scm.duration
                fs = scm.source(:, :, which_source);
            else
                fs = zeros(scm.grid_size);
                if ~isnan(scm.post_ictal_source_drive)
                    fs = fs + scm.post_ictal_source_drive;
                end
            end
        end
        function plot_sigmoids(scm, x)
            if nargin < 2, x = linspace(0, 2, 100); end
            figure; ax = axes(figure, 'nextplot', 'add');
            for ff = {'sigmoid_kdVe', 'sigmoid_kdVi', 'sigmoid_kD'}
                y = scm.sigmoid(x, scm.(ff{:}));
                plot(ax, x, y, 'displayname', ff{1}(9:end))
            end
            legend(ax)
        end
        function clean(scm)
            delete(sprintf('%s_%d_*.mat', scm.basename, scm.sim_num));
        end
        function map_ = generate_contagion(scm)
            % pre-generate the contagion style recruitment map
            rng(scm.sim_num+1);  % add 1 in case you use sim_num=0
            
            map_ = inf(scm.grid_size);
            if scm.expansion_rate <= 0, return; end
            dt_ = scm.dt * 10;  % lower precision is ok; worth increase in speed
            p_wavefront = 2^(scm.expansion_rate * dt_/scm.dx) - 1;  % area is recruited at this rate
            
            recruited = false(scm.grid_size);  % set state negative to indicate "not recruited"
            recruited(scm.stim_center(1), scm.stim_center(2)) = true;  % stim center is recruited at t = 0;
            map_(recruited) = 0;
            tt = dt_;  % advance tt
            
            % repeat until end of simulation or all nodes are recruited
            while tt < (scm.duration + scm.padding(2)) && any(isinf(map_), 'all')
                                
                % find non-recruited points with recruited neighbors
                boundary = conv2(recruited, [0 1 0; 1 -4 1; 0 1 0], 'same') > 0;  
                
                % advance wavefront randomly
                dice = rand(size(recruited));  
                wavefront = (dice < p_wavefront) & boundary;
                
                recruited(wavefront) = true;
                
                map_(wavefront) = tt;
                
                tt = tt + dt_;
            end

        end
        function map_ = generate_map(scm)
            % Generates a linear radial distance from source map with some
            % 2D gaussian noise
            map_ = inf(scm.grid_size);
            if scm.expansion_rate <= 0, return; end
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
        function dims = get.dimsNP(scm)
            dims = arrayfun(@min, scm.dimsNP, scm.grid_size);
        end
        function ctr = get.centerNP(scm)
            if isempty(scm.centerNP)
				scm.centerNP = scm.grid_size / 2; 
            end
            ctr = scm.centerNP;
            xx = 1:scm.grid_size(1);
            yy = 1:scm.grid_size(2);
            [~, ext_x] = sort(abs(xx - ctr(1))); ext_x = sort(xx(ext_x(1:scm.dimsNP(1))));
            [~, ext_y] = sort(abs(yy - ctr(2))); ext_y = sort(yy(ext_y(1:scm.dimsNP(2))));
            ctr = mean([ext_x', ext_y']);

        end
        function [inds, ext_x, ext_y] = NPinds(scm)
            % Computes the indices extracted as the MEA. Recomputes
            % scm.centerNP and scm.dimsNP when called.
            
            mea_mask = false(scm.grid_size);
            ctr = scm.centerNP;
            
            xx = 1:scm.grid_size(1);
            yy = 1:scm.grid_size(2);
            [~, ext_x] = sort(abs(xx - ctr(1))); ext_x = sort(xx(ext_x(1:scm.dimsNP(1))));
            [~, ext_y] = sort(abs(yy - ctr(2))); ext_y = sort(yy(ext_y(1:scm.dimsNP(2))));
            mea_mask(ext_x, ext_y) = true;
            inds = find(mea_mask);
    
            
        end
        function inds = ECinds(scm)
            % Never used this... might need help
            x_offset = ( (1:scm.dimsEC(1)) - floor(scm.dimsEC(1)/2) ) * scm.scaleEC;
            y_offset = ( (1:scm.dimsEC(2)) - floor(scm.dimsEC(2)/2) ) * scm.scaleEC;
            [xx, yy] = ndgrid(x_offset, y_offset);

            inds = sub2ind(scm.grid_size, ...
                scm.centerNP(1) + xx(:), ...
                scm.centerNP(2) + yy(:));
        end
        function h = show_layout(scm)
            
            X = zeros(scm.grid_size);
            X(logical(scm.IZ)) = 1;
            X(scm.NPinds) = 4;  % mea
            X(logical(scm.source)) = 2;
            
            X(scm.stim_center(1), scm.stim_center(2)) = 3;  % IW
%             [xx, yy] = ind2sub(scm.grid_size, scm.NPinds);
%             inds = sub2ind(scm.grid_size, yy, xx);
            
            
%             X(~scm.excitability_map) = -1;
            h = figure('name', 'layout'); imagesc(X);
            cb = colorbar;
            set(cb, 'ticks', 1:4, 'limits', [0 4], ...
                'ticklabels', {'IZ', 'FS', 'IWs', 'MEA'})
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
            dx_ = scm.centerNP(1) - scm.stim_center(1);
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
    
    properties 
        Qe_movie
        Ve_movie
        K_movie
    end
	
	properties  % updates
        
        
        % Epileptogenic zone and fixed source. A baseline increase in excitability over a
        % subregion of the sim. Default shapes set in the getters
        IZ = 1.2
        
        % A predetermined IW map. Generate this if the IW is independent of
        % the rest of the sim. It should be a grid of recruitment times. I
        % tried using this with a t:=Cr function (with 2D gaussian noise;
        % see <SCM.generate_map>) and the recruitment pattern looked pretty
        % good, but the smoothness of it made it so no TW formed after the
        % IW until the FS was inside the recruited region
        map  
        
        % Parameters for the inhibitory collapse portion of the IW. The
        % model assumes an increase in inhibition (depression of dVe)
        % followed by inhibitory collapse, modeled here as a drop in the
        % max inhibitory firing rate. 
        %     Used in: SCM.Qi_max_fun
        %     Elements: 
        %       drop: the amount by which the firing rate decreases
        %       offset: time (s) of peak drop relative to duration of dVe depression
        %       width: sd of gaussian used to model the drop 
        % 
        Qi_collapse = [20 1 .5]  % [drop, offset, width]
        
        
        % SIGMOIDS
        % Change dV* and Dii from dynamic functions of [K+]o to direction
        % sigmoid response functions. The dynamic functions had no
        % decrease so dV*/Dii were changing exponentially with K+. I like
        % the idea of a direct (logistic function sigmoid) relationship
        % better. I think the hypothesis is the K+ drives changes in
        % voltage (and gap junction/ephaptic coupling) directly, rather
        % than changing the rate of change of these functions? 
        %
        % Don't use xxx_center/xxx_width style anymore. Use sigmoid_xxx instead
        % i.e. sigmoid_xxx = [xxx_center params.xxx(2) xxx_width*4]
        
% 		kdVe_center = 0.8  % center of K-->dVe sigmoid  [Martinet = 0.8]
% 		kdVe_width = .1  % ... and width  [Martinet = 0.8]
% 		kD_center = 3e-4  % center of K-->Dii sigmoid [Martinet = .06 / .85]
% 		kD_width =  5e-4  % ... and width [Martinet = .01 / .06]
        
        gap_resting_update {mustBeMember(gap_resting_update,{'sigmoid','dynamic'})} = 'sigmoid'
        sigmoid_kdVe = [0.5 .7 10]  % [mid max slope] (slope is ~max/width*4 if you prefer to think of it that way)
        sigmoid_kdVi = [.5 .3 10]
        sigmoid_kD = [1 .3 -4]
        
        
        % Defines the amount of time the IW spends on each electrode.
        % Specifically, defines the length of time (in seconds) that there
        % is a change to dVe during the IW. Additionally, the timing of the
        % inhibitory collapse is based on this as well (with an additional
        % offset parameter)
        excitability_map  % boundary to IW spread
        
        
        % Create an inhibitory IW. 
        drive_style (1,:) char = 'inhibitory'  % Allow inhibitory driving (EDS, 1/26/21)
        I_drive (1, 1) double = .7  % Allow inhibitory driving (EDS, 1/26/21)
        post_ictal_I_drive (1, 1) double = nan  % Allow inhibitory driving (EDS, 1/26/21)
        
        
        % For simplicity, remove (for now) the depolarization block
        % function that was added in Martinet. With this parameter set, the
        % voltages never get high enough for it to take effect. Maybe look
        % into this more later.
        depo_block = 0  % apply depolarization block
        
    end
    
    properties  % original-ish
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
		subsample (1,1) double {mustBePositive} = Inf  % Allow downsampling when creating mea
		return_fields (1,:) = {'Qe', 'Ve'}  % Qe required to make mea
		out_vars (1,:) = {'Qe', 'Ve'}  % Define which variables you would like to visualize (can be any in IC)
		source % Define fixed, alternating sources of excitation
% 		seed = rng  % Set seed for repeatability
		label (1,:) char = 'SCM'
		t0  % used to keep track of progress in sims 
        
        % model
		stim_center (1,2) = [20 20]
		grid_size (1,2) {mustBePositive} = [50 50] % size of grid to simulate (must be even)
		dt = 2e-4
		dx = .4  % (cm) 
		expansion_rate (1,1) double {mustBeNonnegative} = .625  % in cm^2/s; set to 0 for fixed source
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
		K (1,2) double = [-Inf Inf] % extracellular potassium concentration (cm^2)
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
    
    
    methods  % Emily parameter functions
        function FS = get.source(scm)
            % Default the fixed source to a .25x.25 cm ellipse centered at
            % [2.2, 2.2] cm with the value given in scm.FS (this is in
            % addition to IZ)
            if isempty(scm.source), scm.source = 0; end
            if numel(scm.source) == 1
                center = round( [2.2 2.2] / scm.dx );
                radius = round( [.25 .25] / scm.dx );
                scm.source = double(scm.ellipse(center, radius)) * scm.source;
            end
            FS = scm.source;
        end
        function IZ = get.IZ(scm)
            % Default the epileptogenic zone to an ellipse with the value
            % given in scm.IZ
            if numel(scm.IZ) == 1
                scm.IZ = double(scm.ellipse()) * scm.IZ;
            end
            IZ = scm.IZ;
        end
        function qm = Qi_max_fun(scm, state)
            % Model of inhibitory collapse. Max inhibitory firing rate
            % collapses following an inverse gaussian (this was chosen
            % only because it gives a smooth drop and recovery; no proposed
            % relation to mechanism). 
            if nargin < 2 || all(isnan(scm.Qi_collapse)) % if no state is given, return Qi_max
                qm = scm.Qi_max;
            else  % else update Qi_max relative to IW state
                em = scm.excitability_map;
                em(em <= 0) = -inf;
                qm = scm.Qi_max ...
                    - scm.Qi_collapse(1) * scm.gaussian( ...
                        state, em + scm.Qi_collapse(2), scm.Qi_collapse(3));
            end 
        end
        function map = get.excitability_map(P)
			if isempty(P.excitability_map)
                map = zeros(P.grid_size);
                map(P.ellipse) = .5;
                P.excitability_map = map;
			end
            map = P.excitability_map;

			assert(all(size(map) == P.grid_size))
		end

    end
    
	methods  % original model parameter getters
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
    
    
    methods  % functions for convenience/setup
        function str = mea_path(scm)
            str = sprintf('%s/%s/%s_Seizure%d_Neuroport_%d_%d.mat', ...
                pwd, scm.label, scm.label, scm.sim_num, scm.padding);
        end
        function map = ellipse(scm, C, R, value)
            % map = scm.ellipse(center=grid_size/2, radius=grid_size/2-4)
            % Returns a binary map matching scm.grid_size with an ellipse
            % centered at C (1x2) with radius R (1x2)
            
            if nargin < 2 || isempty(C),  C = scm.grid_size ./ 2; end
            if nargin < 3 || isempty(R); R = floor(scm.grid_size / 2 - 4); end
            if nargin < 4 || isempty(value); value = 1; end
            
            if numel(R) == 1, R = [R R]; end
            [ix, iy] = ind2sub(scm.grid_size, 1:prod(scm.grid_size));

            map = zeros(scm.grid_size);
            ellipse = sum( ([ix' iy'] - C).^2 ./ R.^2, 2 ) < 1;
            map(ellipse) = value;
        end
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

		function dt = get.dt(P)
			if P.dt > 2e-4 && (any(P.IC.Dii(:) >= 0.87) || P.SS.v(:) > 140)
				P.dt = 2e-4;
				warning('High Dii or SS.v; setting dt to 2e-4.');
			end
			dt = P.dt;

        end
        function outs = get.out_vars(scm)
            outs = [scm.out_vars, {'Qe', 'Ve'}];
            [~, o_inds] = unique(outs);
            outs = outs(sort(o_inds));
            
        end
		function grid_size = get.grid_size(p)
			if any(mod(p.grid_size, 2)), p.grid_size = p.grid_size + mod(p.grid_size, 2); end
			grid_size = p.grid_size;
        end
%         function dims = get.dimsNP(scm)
%             dims = arrayfun(@min, scm.dimsNP, scm.grid_size);
%         end
        
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
		function y = sigmoid(x, p)
            y = p(2) ./ (1 + exp(-p(3) * (x - p(1))));
        end

        function y = gaussian(x, mu, sigma)
            y = exp(-(x - mu).^2 ./ (2 .* sigma.^2));
        end
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

