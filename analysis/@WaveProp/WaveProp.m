classdef WaveProp
	
	properties 
		Name
        Patient
        Seizure
		Vx = nan
		Vy = nan
		p = nan
		Beta = nan(1, 6)
		Curvature = nan
		Data
		sig = 0.05
		NClust
		ClustSize
		Quadratic = false
		RotateBy = 0
        MinFinite = 30
% 		HalfWin
% 		FBand
    end
    
    
    properties (Transient = true, Hidden = false)
        Original = false
    end
	
	properties (Transient = true, Hidden = true)
% 		t_inds
% 		res
		
		Magnitude
        AltMagnitude
		Direction
		Time
		Position
		scale_quiver = 1
		complexZ
		Inds = []
		Early = [-Inf 25]
		Late = [30 Inf]
		MinDetections = 100
	end
	
	properties (Hidden = true)
		WPParamNames = "Quadratic"
		t0
	end
	
	properties (Dependent = true, Hidden = false)

		Z
        time
		mask
		logp
		NormedMagnitude
% 		Phi_mean_early
% 		Phi_mean_late
		Phi_std_early
		Phi_std_late
		First_detection
		N_detections_early
		N_detections_late
        FName
    end
    
    methods  % imports
        [H, centers, T] = hist(obj, window, mean_center, show_directions, h)
        ax = particle_smoother(W, ax, NPARTS, SPREAD_FACTOR, SD_WIN, ANG_RES, SMOOTHING_WIN, SHOW_DIRECTIONS, SHOW_EVOLUTION)
    end
	
	methods  % Getters and Setters
        
        function fname = get.FName(self)
            if self.Original
                fname = sprintf('wave_prop/%s_Seizure%s_Neuroport_10_10_wave_prop.mat', self.Patient, self.Seizure);
            else
                fname = sprintf('%s_Seizure%s_fits.mat', self.Patient, self.Seizure);
            end
        end
        
        function pat = get.Patient(self)
            if isempty(self.Patient)
                info = strsplit(self.Name, {'_', 'Seizure'});
                self.Patient = info{1};
                self.Seizure = info{2};
            end
            pat = self.Patient;
        end
        
        function sz = get.Seizure(self)
            if isempty(self.Seizure)
                info = strsplit(self.Name, {'_', 'Seizure'});
                self.Patient = info{1};
                self.Seizure = info{2};
            end
            sz = self.Seizure;
        end
        
		function inds = get.Inds(s)
			inds = s.Inds;
			if isempty(inds)
				inds = (1:numel(s.t0))';
			end
        end
        
        function self = set.time(self, value)
            self.t0 = value;
        end
        
		function obj = set.Beta(obj, val)
			assert(isempty(val) || any(size(val) == 6));
			if isempty(val) || size(val, 2) == 6
				obj.Beta = val;
			else
				obj.Beta = val';
			end
        end
        
		function beta = get.Beta(obj)
			beta = obj.Beta(obj.Inds, :);
        end
        
		function mgn = get.NormedMagnitude(obj)
			mgn = obj.Magnitude;
			mgn(isoutlier(mgn)) = nan;
			mgn = (mgn - nanmean(mgn)) / nanstd(mgn);
			mgn = mgn(obj.Inds);
        end
        
		function logp = get.logp(s)
			logp = -log(s.p)/log(10);
			logp = logp(s.Inds);
        end
        
		function Z = get.complexZ(s)
			Z = complex(s.Vx, s.Vy);
			Z = Z(s.Inds);
        end
        
        function minfinite = get.MinFinite(s)
            if numel(s.MinFinite) > 1
                s.MinFinite = max(s.MinFinite);
            end
            minfinite = s.MinFinite;
        end
        
		function mask = get.mask(s)
			M = sqrt(s.Vx.^2 + s.Vy.^2);
			M(M > 1e6) = nan;
			M(isoutlier(M)) = nan;
			M = M(s.Inds);
            num_finite = sum(isfinite(s.Data), [2 3]);
            M(num_finite < s.MinFinite(1)) = nan;
			mask = s.p > s.sig | isnan(M);
            
            % mask times when VNS is active in CUCX5
            if strcmpi(s.Name, 'cucx5_seizure3')
                vns_times = s.time >= 13 & s.time <= 29;
                mask = mask | vns_times;
            end
            if strcmpi(s.Name, 'cucx5_seizure6')
                vns_times = s.time >= 46 & s.time <= 62;
                mask = mask | vns_times;
            end
        end
        
		function D = get.Data(s)
			if size(s.Data, 1) ~= numel(s.t0)
				time_dim = find(size(s.Data) == numel(s.t0));
				if isempty(time_dim)
					D = reshape(s.Data, [numel(s.t0) size(s.Data)]);
				else
					D = shiftdim(s.Data, time_dim-1);
				end
			else
				D = s.Data;
			end
			D = D(s.Inds, :, :);
        end
        
		function time = get.time(s)
			time = s.t0;
			time = time(s.Inds);
        end
        
		function Z = get.Z(s)
			Z = s.Direction;
			Z = Z(s.Inds);
        end
        
        function C = get.Curvature(s)
            if ~s.Quadratic, C = 0; return; end
			C = s.Curvature;
            if numel(C) == 1 && numel(s.Inds) > 1
                C = C * ones(size(s.Inds));
            end
			C(s.mask) = nan;
			C = C(s.Inds);
        end
        
		function RB = get.RotateBy(s)
			RB = s.RotateBy .* ones(size(s.Vx));
			RB = RB(s.Inds);
        end
        
		function D = get.Direction(s)
			D = atan2(s.Vy, s.Vx);
			D = D(s.Inds);
			D(s.mask) = nan;
			D = D - s.RotateBy;
			D = angle(exp(1j * D));  % Keep in range [-pi pi]
        end
        
		function M = get.Magnitude(self)
            % in mm/s
			M = sqrt(self.Vx.^2 +  self.Vy.^2) * .4;  % Fits return electrodes/second. Convert to mm/s
			M = M(self.Inds);
            M(self.mask) = nan;
        end
        function M = get.AltMagnitude(self)
            t_ext = range(self.Data, [2 3]);  % time taken for wave to cross MEA
            M = 5.67 ./ t_ext;  % speed in mm/second (5.67 mm = length of MEA diagonal)
            M = M(self.Inds);
            M(self.mask) = nan;
        end
        
		function psig = get.p(s)
			psig = s.p(s.Inds);
		end
		
		%% Summary stats
		function e = early(s)
			e = s.time > s.Early(1) & s.time < s.Early(2);
		end
		function l = late(s)
			l = s.time > s.Late(1) & s.time < s.Late(2);
		end
		function [phimn, conf] = Phi_mean_early(s)
            if s.N_detections_early < s.MinDetections, [phimn, conf] = deal(nan); return; end
			[phimn, ul] = circ_mean(s.Direction(s.early), [], [], 'omitnan');
			conf = ul - phimn;
			phimn = angle(exp(1j*phimn));  % keep in range [-pi pi]
		end
		
		function [phimn, conf] = Phi_mean_late(s)
            if s.N_detections_late < s.MinDetections, [phimn, conf] = deal(nan); return; end
			[phimn, ul] = circ_mean(s.Direction(s.late), [], [], 'omitnan');
			conf = ul - phimn;
			phimn = angle(exp(1j*phimn));  % keep in range [-pi pi]
		end
		function phisd = get.Phi_std_early(s)
			phisd = circ_std(s.Direction(s.early), [], [], 'omitnan');
			phisd(s.N_detections_early < s.MinDetections) = nan;
		end
		function phisd = get.Phi_std_late(s)
			phisd = circ_std(s.Direction(s.late), [], [], 'omitnan');
			phisd(s.N_detections_late < s.MinDetections) = nan;
		end
		function t1 = get.First_detection(s)
			times = [s.time; nan];
			% t1 = times(find(isfinite(s.Direction), 1));
			t1 = quantile(times(isfinite(s.Direction)), .025);
		end
		function nd_early = get.N_detections_early(s)
			nd_early = numfinite_(s.Direction(s.early));
		end
		function nd_late = get.N_detections_late(s)
			nd_late = numfinite_(s.Direction(s.late));
		end
		
	end
	
	methods 
		
		%% Plotting
		function ax = plot(obj, varargin)
			% ax = plot(obj, type='2D', t0, ax); % if object has only one time point, t0=obj.time
			directive = '';
			chars = cellfun(@ischar, varargin);
			if any(chars), directive = varargin{chars}; varargin(chars) = []; end
			switch upper(directive)
				case '3D'
					ax = obj.plot3D(varargin{:});
				otherwise
					ax = obj.plot2D(varargin{:});
			end
        end
        
		function ax = plot3D(obj, varargin)
			% ax = plot(obj, t0, ax); % if object has only one time point, t0=obj.time
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
			if numel(obj.time) == 1
				F = obj;
				units = '';
			else				
				[~, ind] = min((obj.time - tt).^2);
				F = obj.sub(ind);
				F.Magnitude = obj.NormedMagnitude(ind);
				units = 'std';
			end
			data = squeeze(F.Data);
			N = size(data);
			[p1, p2] = ind2sub(N, find(isfinite(data)));
			scatter3(ax, p1, p2, F.Data(isfinite(data)), [], F.Data(isfinite(data)), 'filled'); 
			hold(ax, 'on');
			[xx, yy] = ndgrid(1:N(1), 1:N(2));
			zfit = F.Beta(1) + F.Beta(2) * xx + F.Beta(3) * yy;
			im = surf(ax, xx, yy, min(F.Data(:)) * ones(N), data, 'linestyle', 'none'); 
			im.Tag = 'Data';
			surfax = surf(ax, xx, yy, zfit, 'linestyle', 'none', 'facealpha', .5);
			surfax.Tag = 'Fit';
			title(ax, sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', F.time, F.p, F.Magnitude, units)); 
			hold(ax, 'off');
			ax.Tag = checkname(['figs' filesep obj.Name '_' class(obj) '_' num2str(F.time) '_plot3D']);
			
        end
        
		function ax = plot2D(s, varargin)
			% ax = plot2D(obj, tt, ax)
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
			if ndims(s.Data) == 3
				[~, which_t] = min(abs(tt - s.time));
				data = squeeze(s.Data(which_t, :, :));
				vx = s.Vx(which_t);
				vy = s.Vy(which_t);
			else
				data = s.Data;
				vx = s.Vx;
				vy = s.Vy;
			end
			imagesc(ax, data');
			axis(ax, 'xy'); 
			set(ax, 'xtick', [], 'ytick', [])
			hold(ax, 'on')
			X = (size(data) + 1) / 2;
			quiver(ax, X(1), X(2), vx, vy, s.scale_quiver, ...
				'LineWidth', 2, ...
				'MaxHeadSize', 1, 'color', [1 0 0])
			hold(ax, 'off')
			ax.Tag = checkname(['figs' filesep s.Name '_' class(s) '_' num2str(tt) '_plot2D']);
        end
        
		function [d, xi, bw] = ksdensity(obj, ax, rotateby, xres, bw, varargin)
			% [d, xi, bw] = ksdensity(ax=gca, rotateby=0, xres=pi/500, bw=.15*pi/3, ::plot directives::)
			if nargin < 5 || isempty(bw), bw = .15*pi/3; end
			if nargin < 4 || isempty(xres), xres = 500; end
			if nargin < 3 || isempty(rotateby), rotateby = 0; end
			if nargin < 2 || isempty(ax), ax = gca; end
			if rotateby ~= 0
				temp = angle(exp(1j * (obj.Direction - rotateby)));
			else
				temp = obj.Direction;
			end
			temp = [temp(:) - 2*pi; temp(:); temp(:) + 2*pi];
			gridx2 = linspace(-3*pi, 3*pi, 6*xres + 1);
			[d, xi, bw] = ksdensity(temp, gridx2, 'bandwidth', bw);
			plot(ax, xi, d, varargin{:});
			xlim(ax, [-pi pi])
		end
		
		%% Main
		function obj = compile_results(obj, varargin)
			% obj = compile_results(::'nocluster'::)
            
			CLUSTER = true;
			if ismember('nocluster', varargin), CLUSTER = false; end
			temp = nan(max(obj.Position));
			locs = sub2ind(size(temp), obj.Position(:, 1), obj.Position(:, 2));
			if numel(unique(obj.Data(isfinite(obj.Data)))) < 3
				temp(locs) = obj.Data;
				obj.Data = temp;
				return
			end
			warning('off', 'stats:clusterdata:MissingDataRemoved')
			if CLUSTER
				data = obj.use_largest_cluster(obj.Data, obj.Position);
			else 
				data = obj.Data;
			end
			[V, obj.p, beta] = obj.fit_data(data, obj.Position, obj.Quadratic);
			obj.Vx = V(1);
			obj.Vy = V(2);
			obj.Beta = beta;		
			obj.Curvature = WaveProp.curvature_from_beta(beta);
			ZZ = complex(V(1), V(2));
			ZZ(abs(ZZ) > 1e6) = nan;
			obj.Magnitude = abs(ZZ);
			obj.Direction = angle(ZZ);
			temp(locs) = data;
			obj.Data = single(temp);
        end
        
		function obj = parse_inputs(obj, varargin)
			for ii = 1:2:numel(varargin)
				ff = validatestring(varargin{ii}, [obj.ParamNames(:); obj.WPParamNames(:)]);
				obj.(ff) = varargin{ii+1};
			end
        end
        
        function cst = consistency(obj, W)
            % cst = obj.consistency(W=1000);
            % Default to 1000 ms window
            if nargin < 2, W = 1000; end  
            k = round(W/diff(obj.time(1:2)) / 1e3);  % convert window from ms to samples
            c = cos(obj.Direction);
            s = sin(obj.Direction);
            cS = movmean(c, k, 'omitnan');
            sS = movmean(s, k, 'omitnan');
            
            cst = vecnorm([cS sS], 2, 2);
            N = movsum(isfinite(c), k);
            cst(N < 5) = nan;
            
        end
        
        
		
		%% Updates
		function sN = refit_data(s)
			W = warning;
			warning('off', 'stats:statrobustfit:IterationLimit');
			s.Inds = [];  % Remove time resampling
			data = s.Data;
			eval(sprintf('sN = %s();', class(s)));
			sN = s;
			sN.Beta = [];
			locs = find(isfinite(data));
			[p1, p2] = ind2sub(size(s.Data), locs);
			sN.Data = data(locs)';
			sN.Position = [p1(:) p2(:)];
			sN = sN.compile_results('nocluster');
			warning(W);
        end
        
		function s = compute_curvature(s)
			S = warning;
			warning('off', 'stats:statrobustfit:IterationLimit');
			INDS = s.Inds; s.Inds = [];  % temporarily remove time filter
			
			N = numel(s.time);
			s.Beta = nan(N, 6);
			s.Curvature = nan(N, 1);
			
			dims = size(s.Data);
			[xx, yy] = ndgrid(1:dims(2), 1:dims(3));
			position = [xx(:) yy(:)];
			for ii = 1:N
				data = squeeze(s.Data(ii, :, :));
				finite = isfinite(data(:));
				[V, p0, beta] = WaveProp.fit_data(data(finite), position(finite, :)); 
				s.Vx(ii) = V(1); 
				s.Vy(ii) = V(2); 
				s.p(ii) = p0; 
				s.Beta(ii, :) = beta; 
% 				Zc = complex(V(1), V(2));
% 				s.Magnitude(ii) = abs(Zc);
% 				s.Direction(ii) = angle(Zc);
				s.Curvature(ii) = WaveProp.curvature_from_beta(beta);
			end
			warning(S);
			s.Inds = INDS;  % re-apply time filter
        end
        
		function [mn, ci] = mean(s, data)
			if nargin < 2, data = s.Direction; end
            if all(isnan(data)), mn = nan; ci = nan; return; end
			[mn, ul] = circ_mean(data, [], [], 'omitnan');
			ci = ul - mn;
        end
        
        function sd = std(s, data)
            if nargin < 2, data = s.Direction; end
            if all(isnan(data)), sd = nan; return; end
            sd = circ_std(data(~isnan(data)));
        end
        
		function F = sub(obj, ind)
			obj.Inds = [];
			N = min(numel(obj.time), numel(obj.Vx));
			eval(sprintf('F = %s;', class(obj)));
			F.t0 = obj.t0(ind);
			for f = string(fieldnames(obj)')
				temp = obj.(f);
				sz = size(temp);

                try
				if all(sz < N) || ischar(obj.(f))
					F.(f) = temp;
				elseif isvector(temp)
					F.(f) = temp(ind);
				elseif ismatrix(temp)
					F.(f) = temp(ind, :);
				else
					assert(size(temp, 1) >= N);
% 					dim = 3;
% 					temp = shiftdim(temp, dim-1);
					sz = size(temp);
					temp = reshape(temp, N, []);
					data = reshape(temp(ind, :), [1 sz(2:end)]);
					F.(f) = data;
                end
                catch ME
                    if ~strcmpi(ME.identifier, 'MATLAB:class:noSetMethod')
                        rethrow(ME)
                    else
                        continue
                    end
                end

			end
        end
        
		function F = split(obj)
			N = min(numel(obj.time), numel(obj.Vx));
			F(N) = obj.sub(N);
			for ii = 1:N
				F(ii) = obj.sub(ii);
			end
        end
        
		function d = diff(s, other)
			other = other.resample_t0(s.time);
% 			t_inds = interp1(other.time, 1:length(other.time), s.time, 'nearest', 'extrap');
			d = angle(exp(1j*(s.Direction - other.Direction)));
		end
		
		%% Reshaping
		function [obj, iq] = resample_t0(obj, t_new)
			iq = interp1(obj.t0, 1:numel(obj.t0), t_new, 'nearest', 'extrap');
			obj.Inds = iq(:);
		end	
		function S = compile(obj, t_fun)
			% Compile an array of WaveProp objects into one
			% S = obj.compile(t_fun=@(t) t0 < Inf)
			if nargin < 2, t_fun =@(t0) t0 < Inf; end
			S = obj(1);
			for ii = 1:numel(obj)
				t_mask = t_fun(obj(ii).time);
				obj(ii).p(~t_mask) = Inf;
			end
% 			S.p = cat(1, obj.p);
% 			S.Magnitude = cat(1, obj.Magnitude);
			for ff = string(fieldnames(obj)')
				
% 				if ismember(lower(ff), {'name'}), continue, end
% 				if strcmpi(ff, 'data'), S.(ff) = cat(3, obj.(ff)); continue; end
				try
					S.(ff) = cat(1, obj.(ff));
				catch ME
					S.(ff) = {obj.(ff)};
					if ~strcmpi(ME.identifier, 'MATLAB:catenate:dimensionMismatch')
						rethrow(ME)
					end
				end
			end
% 			t_mask = t_fun(S.time);
% 			S.p(t_mask) = Inf;
			t0d = [0; diff(S.time)];
			t0d(t0d <= 0) = 1;
			S.time = cumsum(t0d);
			S.sig = min(S.sig);
			
		end
		
	end
		
	methods (Static)
		function [ax, tt, directive] = parse_plot_inputs(varargin)
			ax = []; tt = 0;
			for arg = varargin
				if isa(arg{:}, 'matlab.graphics.axis.Axes')
					ax = arg{:};
				elseif isnumeric(arg{:})
					tt = arg{:};
				elseif ischar(arg{:})
					directive = arg{:};
				else
					error("Unrecognized input")
				end
            end
            if isempty(ax), ax = axes(figure); end
		end
		
		function [mat, times] = WP2mat(obj, field, times)
			% [mat, tt] = WP2mat(S, field, times=[])
			% If no times are given, the time points from the first metric
			% in obj are used
			if nargin < 3, times = []; end
			X = {};
			for ff = string(fieldnames(obj)')
				if isempty(times)
					times = obj.(ff).time; 
				else
					obj.(ff) = obj.(ff).resample_t0(times);
				end
% 				[times, ii] = unique(S.(ff).time);
% 				iq = interp1(times, ii, tq, 'nearest', 'extrap');
				X = [X {obj.(ff).(field)}]; %#ok<AGROW>
			end
			mat = cat(2, X{:});
		end
		function fits = file2fit(res)
			
			f = string(fieldnames(res));
			for ff = f'
				m = rename_metrics(ff);
				fits.(m) = WaveProp.struct2obj(res.(ff));
			end
			
		end
		
		function res = cell2mat(data, position)
			N = numel(data);
			assert(N == numel(position));
			dims = max(cell2mat(position));
			res = nan(dims(1), dims(2), N);
			for ii = 1:N
				temp = nan(dims);
				P = position{ii}; D = data{ii}(:);
				if isempty(P), continue, end
				locs = sub2ind(dims, P(:, 1), P(:, 2));
				[G, l] = findgroups(locs);
				value = splitapply(@mean, D, G);
				temp(l) = value;
				res(:, :, ii) = temp;
			end
			
		end
		
		function C = curvature_from_beta(beta)
% 			A = sum(beta(4:5));
% 			B = sum(beta([4 6]));
			A = beta(5);
			B = beta(6);
			C = complex(A, B);
		end
		
		function [b, P0 ] = estimate_wave( data, position, varargin )
			%ESTIMATE_WAVE Attempts to fit a two-dimensional plane to the delays between electrodes
			%organized in space based on their positions.
			%   [SRC_DIR,SPEED,CI_DIR,CI_SP]=ESTIMATE_WAVE(DELAY,POSITION,VARARGIN)
			%   fits a plane to the DELAYs between all electrodes and the most central
			%   electrode based on the electrode POSITIONs in space. If there are
			%   enough electrodes with defined delays and the fit is significant, then
			%   the direction SRC_DIR toward the source of the wave propagation is 
			%   returned as well as the SPEED of the wave and confidence intervals
			%   around those values (CI_DIR, CI_SP) computed via bootstrapping.
			%   Optional: if the last input is 'plot', a representation of the delays
			%   at each electrode position with the fitted plane will be shown.

			MIN_RATIO_FINITE = 0.05;   % we require 5% of electrodes to have a defined delay (this is overwritten later by MinFinite)

			P0 = nan;
			b = nan;
			
			finite = numfinite_(data);
			
			
			if finite > 3 && finite/numel(data) >= MIN_RATIO_FINITE  % check enough delay data is not NaN.
				if any(ismember({'quadratic', 'quad'}, lower(varargin)))
					X = [position prod(position, 2) position.^2];  % predictors
				else 
					X = position;
				end
				[b,stats] = robustfit(X, data, 'fair');                     % fit the delay vs two-dimensional positions
% 				H = [0 1 0; 0 0 1];  % These are default values
% 				c = [0 ; 0];
				P0 = linhyptest(b, stats.covb, [], [], stats.dfe);  % perform F test that last two coefs are both 0.
            end
            
        end
        
		function [V, p, beta] = fit_data(data, position, quadratic)
			if nargin < 3, quadratic = false; end
			V = nan(1, 2);
			beta = nan(1, 6); 
			[b, p] = WaveProp.estimate_wave(data, position);
			[betaQ, pQ] = deal(nan);
			if quadratic
				try
					[betaQ, pQ] = WaveProp.estimate_wave(data, position, 'quad');
				catch ME
					if ~strcmpi(ME.identifier, 'stats:robustfit:NotEnoughData')
						rethrow(ME)
					end
					pQ = nan;
				end
			end
			if isnan(p) && isnan(pQ), return; end
			if isnan(pQ), beta(1:3) = b; beta(4:end) = 0;
			elseif isnan(p), beta = betaQ; p = pQ; 
			elseif pQ < p; beta = betaQ;
			else, beta(1:3) = b; beta(4:end) = 0;			
			end
% 			if numel(beta) == 3, beta = [beta; nan(3, 1)]; end
			% beta is invalid if nan or if the slope is 0 in both directions
			invalid = all(beta(2:3).^2 < eps);
			if invalid
				return
			end
% 			beta = circshift(beta, -1);
			V = pinv(beta(2:3));
        end
        
		function [reduced_data, T] = use_largest_cluster(data, position)
			% Use the largest cluster of data points for lfp methods
			warning('off', 'stats:clusterdata:MissingDataRemoved');
			T = clusterdata([data(:) position], ...
				'distance', 'seuclidean', ... % standardized euclidean distance
				'cutoff', 1, 'criterion', 'distance');  % Separate large gaps (this distance is somewhat arbitrary right now...)
			[~, largest] = max(histcounts(T, max(T)));
			mask = T == largest;
			reduced_data = data;
			reduced_data(~mask) = nan;
			
% 			dataS = sort(data(:));
% 			if all(isnan(dataS)), reduced_data = data; return; end
% 			dataS(isnan(dataS)) = [];  % excluded nan values
% 			diff_sorted = diff(dataS);  % calculate gaps between nearby data points
% 			big_gaps = isoutlier(unique(diff_sorted));  % find large gaps between datapoints
% 			divides = [0; find(big_gaps); numel(dataS)];  % divisions between clusters in sorted data
% 			cluster_sizes = diff(divides);  % size of clusters
% 			[~, largest] = max(cluster_sizes);  % choose the largest cluster
% 
% 			bounds = dataS(divides(largest+[0 1]) + [1; 0]);
% 			reduced_data = data;
% 			reduced_data(data < bounds(1) | data > bounds(2)) = nan;
        end
        
		function S = resize_obj(s, to_struct)
			if nargin < 2, to_struct = false;  end
			
			S = s(1);
			if to_struct
				S = struct;
				for p = string(properties(s(1))')
					S.(p) = []; 
				end
			end
			S.p = cat(1, s.p);
			S.t0 = cat(1, s.t0);
			S.Vx = cat(1, s.Vx);
			S.Vy = cat(1, s.Vy);
            S.Data = cat(1, s.Data);
			S.Inds = [];
            
% 			S.Magnitude = cat(1, s.Magnitude);
            mco = metaclass(s).PropertyList;
            FF = {mco.Name};
            FF(cat(1, mco.Dependent) | cat(1, mco.Hidden) | cat(1, mco.Transient)) = [];
			for f = string(FF)
				if ismember(lower(f), ...
						{'quadratic', 'sig', 'halfwin', 'fband', ...
						'name', 'minfreq', 'p', 'magnitude', 'Vx', 'Vy', ...
                        'samplingrate', 'patient', 'seizure'})
					continue
				end
				if isempty(s(1).(f)), continue,
                else
                    temp = cat(1, s.(f));
                    try
%                         temp = cat(1, s.(f));
                        t2 = temp;
                        if ismatrix(temp)
                            t2 = unique(temp, 'rows');
                        elseif isvector(S.(f))
                            t2 = unique(temp);
                        end
                        
                        if size(t2, 1) == 1  % if all rows are the same, just keep one
                            S.(f) = t2;
                        else  % otherwise, preserve size
                            S.(f) = temp;
                        end
					    
                    catch ME  
                        % some of the older saves had different sizes. 
                        % Join these as cells and then reshape by hand
                        switch ME.identifier
                            case 'MATLAB:catenate:dimensionMismatch'
                                S.(f) = {s.(f)};
                            case 'MATLAB:class:noSetMethod'
                                continue
                            otherwise
                                rethrow(ME);
                        end
                        
                    end
                    

				end
            end
			
        end
        
		function obj = struct2obj(fstruct)
			% Usage: obj = struct2obj(res)
			
			obj = WaveProp;
			assert(nargin == 1);
			obj.Name = fstruct.Name;
            nn = strsplit(fstruct.Name, {'_', 'Seizure'});
            obj.Patient = nn{1};
            obj.Seizure = nn{2};
			obj.t0 = fstruct.computeTimes/1e3;

			data = fstruct.data;
			position = fstruct.position;
			if iscell(fstruct.data), data = WaveProp.cell2mat(data, position); end
			obj.Data = data;
% 			obj.Position = position;
			try
				obj.NClust = fstruct.n_clust(:);
				obj.ClustSize = fstruct.clust_size(:);
			catch ME
				if ~strcmpi(ME.identifier, 'MATLAB:nonExistentField')
					rethrow(ME)
				end
			end
			obj.Direction = fstruct.Z(:);
			obj.Magnitude = abs(complex(fstruct.V(1, :)', fstruct.V(2, :)'));
			obj.p = fstruct.p';
			obj.Vx = fstruct.V(1, :)';
			obj.Vy = fstruct.V(2, :)';
			beta = fstruct.beta;
			if size(beta, 1) < 6
				beta = [beta; nan(6 - size(beta, 1), size(beta, 2))]';
			end
			obj.Beta = beta;
			
% 			if ~ismember('curvature', lower(fieldnames(fstruct)))
% 				obj = obj.compute_curvature; 				
% 			end
			
		end

		function args = parse_load(varargin)

			P = inputParser;
			p =@(varargin) addParameter(P, varargin{:});

			% Defaults
% 			metrics = {...
% 				'maxdescent', ...
% 				'events', ...
% 				'delays_T10_fband1_13', ...
% 				'delays_T01_fband1_13'}; 
% 			allMetrics = {...
% 				'maxdescent', ...
% 				'events', ...
% 				'delays_T10_fband1_13', ...
% 				'delays_T01_fband1_13', ...
% 				'delays_T10_fband1_50', ...
% 				'delays_T01_fband1_50'}; 

			p('pat', '*');
			p('seizure', '*');
			p('files', []);
% 			p('metrics', metrics, @(c) all(contains(c, allMetrics)));
			p('metrics', []);
			p('sig', 5e-2);
			p('original', false);

			parse(P, varargin{:});

			args = P.Results;
			if isnumeric(args.seizure), args.seizure = num2str(args.seizure); end
			if ischar(args.metrics), args.metrics = {args.metrics}; end
			
			% Cleaning
			if args.original
				fname =@(p, s) sprintf('wave_prop/%s_Seizure%s_Neuroport_10_10_wave_prop.mat', p, s);
			else
				fname =@(p, s) sprintf('%s_Seizure%s_fits.mat', p, s);
			end
			if isempty(args.files)
				if ischar(args.pat) && ischar(args.seizure)
					if strcmpi([args.pat args.seizure], '**')
                        
                        sz = SeizureInfo;
                        pats = sz.patient;
                        seizures = sz.seizure;

                        
						for ii = numel(pats):-1:1
                            fstruct = dir(fname(pats{ii}, num2str(seizures(ii))));
                            if ~isempty(fstruct)
                                files(ii) = fstruct;
                            end
						end
						
						args.files = files;
                        
                        args.seizures = arrayfun(@(ii) ...
                            [sz.patient{ii} '_Seizure' num2str(sz.seizure(ii))], ...
                            (1:numel(pats)), 'uni', 0);
                        
					else
						args.files = dir(fname(args.pat, args.seizure));
                        args.seizures = [args.pat '_Seizure' args.seizure];
					end
                else
                    error('You need to add this functionality: args.seizures not defined for this');
					ii = 1;
					for p = 1:numel(args.pat)
						for s = 1:numel(args.seizure)
							try 
								args.files{ii} = ...
									dir(fname(args.pat{p}, args.seizure{s}));
								ii = ii + 1;
							catch ME
								if ~strcmpi(ME.identifier, 'MATLAB:matrix:singleSubscriptNumelMismatch')
									rethrow(ME)
								end
							end
						end
					end	
				end

			end

        end
        

        function fits = load(mea, metrics)
            % fits = WaveProp.load(mea=[]);
            % Inputs:
            %    mea - can be (1) an MEA class object, (2) a character
            %          string name of a subfolder of WaveFits (i.e.
            %          'MG49_Seizure43'), or (3) it can be empty in which
            %          case fits from all seizures in SeizureInfo will load
            
            if nargin < 2, metrics = []; end
            if nargin < 1 || isempty(mea)
                sz = SeizureInfo;
                metrics = findsharedmetrics_(sz);
                for ii = numel(sz.patient):-1:1
                    names{ii} = sprintf('%s_Seizure%d', sz.patient{ii}, sz.seizure(ii));
                    fits(ii) = WaveProp.load( names{ii}, metrics);
                end
                F.Name = names;
                F.Original = false(size(names));
                for ff = metrics
                    F.(ff{:}) = cat(1, fits.(ff{:})); 
                end
                fits = F;
                
            elseif ischar(mea)
                fits = foldercontents2struct_(['WaveFits/' mea], metrics);
                fits.Name = {mea};  % convert to cell to match previous version of load
            elseif isa(mea, 'MEA')
                fits = foldercontents2struct_(['WaveFits/' mea.Name], metrics);
                fits.Name = {mea.Name};  % convert to cell to match previous version of load
            end
            
            
        end
        
        function fits = COMPLICATEDload(varargin)
            
            
            if nargin > 0 && isa(varargin{1}, 'MEA')
                s = varargin{1};
                varargin(1) = [];
            end
            
            args = WaveProp.parse_load(varargin{:});
            
            if isa(s, 'MEA')
                fname = ['WaveFits/' s.Name];  % folder name
                [fits, metrics] = foldercontents2struct_(fname);
            else
                files = args.files;
				seizures = args.seizures;
				metrics = args.metrics;
				nF = numel(seizures);
				for ii = nF:-1:1
                    fprintf('%d: loading %s\n', ii, seizures{ii});
                    if args.original
                        res = load([files(ii).folder filesep files(ii).name]); 
                        fit_temp = WaveProp.file2fit(res);
                        for mm = metrics
                            fits(ii).(mm{:}) = fit_temp.(mm{:}); 
                        end
						fits(ii).Original = true;
                        
                    elseif ii == nF && isempty(metrics)
                        path = ['WaveFits/' seizures{ii}];
                        [res, metrics] = foldercontents2struct_(path);
                        
                    else
                        res = foldercontents2struct_(path);
                        
						not_member = ~ismember(metrics, fieldnames(res));
						if any(not_member)
							fits = rmfield(res, metrics(not_member));
							metrics(not_member) = [];
						end
					end
					if ~args.original
						res.Name = strrep(seizures{ii}, '.mat', '');
                        res.Original = args.original;
						fits(ii) = res;
                    end

				end	
				for ff = string(fieldnames(fits)')
%                     if strcmpi(ff, 'fits'), continue; end
					if ischar(fits(1).(ff))
						F.(ff) = {fits.(ff)};
					else
						F.(ff) = cat(1, fits.(ff)); 
					end
				end
				fits = F;
			end	
            
        end
        
		function fits = ZZload(varargin)
			s = [];
			if nargin > 0 && isa(varargin{1}, 'MEA')
				s = varargin{1};
				varargin(1) = [];
			end
			args = WaveProp.parse_load(varargin{:});
			if isa(s, 'MEA')
				if args.original
					fname = sprintf( ...
						'wave_prop/%s_Seizure%s_Neuroport_10_10_wave_prop.mat', ...
						s.patient, s.seizure);
				else
					fname = [s.Name '_fits.mat'];
				end
				if ~isempty(args.metrics), fits = load(fname, args.metrics{:});
				else, fits = load(fname); 
				end
				if args.original, fits = WaveProp.file2fit(fits); end
			else
				files = args.files;
				metrics = args.metrics;
				nF = numel(files);
				for ii = nF:-1:1
                    fprintf('%d: loading %s\n', ii, files(ii).name);
                    if args.original
                        res = load([files(ii).folder filesep files(ii).name]); 
                    elseif ii == nF && isempty(metrics)
						res = load([files(ii).folder filesep files(ii).name]);
						metrics = fieldnames(res); 
					else
						res = load([files(ii).folder filesep files(ii).name], metrics{:});
                        
                        
						not_member = ~ismember(metrics, fieldnames(res));
						if any(not_member)
							fits = rmfield(res, metrics(not_member));
							metrics(not_member) = [];
						end
					end
					if args.original
                        fit_temp = WaveProp.file2fit(res);
                        for mm = metrics
                            fits(ii).(mm{:}) = fit_temp.(mm{:}); 
                        end
						fits(ii).Original = true;
					else
						res.Name = strrep(files(ii).name, '.mat', '');
                        res.Original = args.original;
						fits(ii) = res;
                    end

				end	
				for ff = string(fieldnames(fits)')
%                     if strcmpi(ff, 'fits'), continue; end
					if ischar(fits(1).(ff))
						F.(ff) = {fits.(ff)};
					else
						F.(ff) = cat(1, fits.(ff)); 
					end
				end
				fits = F;
			end						
        end
        
		function [out, S, dir_stats] = summary_stats(F, varargin)
			% [out, S, dir_stats] = WP.summary_stats(F, ::'noplot'::, ::'mea'::)
			% Generates stats for each patient
			PLT = true;
			USE_MEA = false;
			if ismember('noplot', varargin), PLT = false; end
			if ismember('mea', varargin), USE_MEA = true; end
			metrics = fieldnames(F);
			mask = ismember(metrics, 'Name');
			metrics(mask) = [];
			Np = numel(F.(metrics{1}));  % number of patients
			Nm = numel(metrics);
			out = table;
			for ii = Np:-1:1
				fname = strsplit(F.Name{ii}, '_');
				names{ii} = [fname{1} ' ' fname{2}(8:end)];
			end
			
			for pp = Np:-1:1
				for mm = string(metrics(:)')
					S(pp).(mm) = F.(mm)(pp);
				end
			end
			out.Weibull_A = nan(Np, Nm);
			for pp = Np:-1:1
				magnitude = WaveProp.WP2mat(S(pp), 'Magnitude');
				out.Magnitude_median(pp, :) = nanmedian(magnitude);
				out.Magnitude_std(pp, :) = nanstd(magnitude);
				
				[Z, T] = WaveProp.WP2mat(S(pp), 'Direction');
				D = WaveProp.compute_dir_stats(Z, T);
				out.Direction_persistence(pp, :) = D.persistence;
				out.Direction_variability(pp, :) = D.dphi_mean;
				for ff = ["phi_mean_early" "phi_mean_late" "phi_std_early" "phi_std_late" "first_detection"]
					out.(ff)(pp, :) = D.(ff);
				end
				dir_stats{pp} = D;
				
				
				for mm = 1:numel(metrics)
					if numfinite_(magnitude(:, mm)) < 10, continue; end
					pd = fitdist(magnitude(:, mm), 'Weibull');
					out.Weibull_A(pp, mm) = pd.A;
					out.Weibull_B(pp, mm) = pd.B;
					
					
				end
				
				if USE_MEA
					whowho = strsplit(names{pp});
					file = dir([whowho{1} filesep '*Seizure' whowho{2} '_*Neuroport_10_10.mat']);
					mea = MEA([file.folder filesep file.name]);
					active = (mea.Time > 0) & (mea.Time < mea.Time(end) - mea.Padding(2));
					out.PeakRatio(pp) = mea.PeakRatio;
					out.MedianLFP(pp) = median(mea.lfp(active, :), 'all');
					out.PeakFR(pp) = quantile(mea.firing_rate(active, :), .95, 'all');
					out.Duration(pp) = mea.Time(end) - mea.Padding(2); 
					ISI = diff(mea.get_wave_times('lfp'));
					pd = fitdist(ISI(:), 'InverseGaussian');
					out.IGmu(pp) = pd.mu;
					out.IGlbda(pp) = pd.lambda;
				end
			end
			
			if PLT
				h = figure; fullwidth(true);
				t = tiledlayout(h, 'flow');
				for ff = string(out.Properties.VariableNames), 
					ax = nexttile(t); stem(ax, out.(ff)); 
					xticks(ax, []); 
					title(strrep(ff, '_', ' ')); 
				end
				set(ax, 'xtick', 1:33, 'XTickLabel', names, ...
					'xticklabelrotation', 90);
				legend(ax, metrics)
			end
			
        end
        
		function out = compute_dir_stats(phi, times, thresh)
			if nargin < 3, thresh = pi/64; end
			[dphi, phi_uw, phi_std] = deal(nan(size(phi)));
			num_metrics = size(phi, 2);
			for jj = num_metrics:-1:1
				phi_temp = phi(:, jj);  % get Z for given patient and metric

				mask = isfinite(phi_temp);  % isolate non-nan
				phi_uw_temp = unwrap(phi_temp(mask));  % unwrap (isolate and then unwrap is questionable ...)
				phi_temp(mask) = phi_uw_temp;
				phi_uw(:, jj) = phi_temp;
				for tt = 1:numel(times)
					mask_tt = (times >= times(tt) - .5) & (times <= times(tt) + .5);
					phi_std(tt, jj) = nanstd(phi_temp(mask_tt));
				end
		% 		phi_temp_sm = smoothdata(phi_temp, 1, 'gaussian', 1, 'SamplePoints', t_temp);  % smooth with a 1 second (gaussian) kernel
		% 		dphi{ii}(:, jj) = [nan; smoothdata(diff(phi_temp_sm))];
				dphi_temp = [nan; smoothdata(diff(phi_temp), 1, 'movmean', 1, 'samplepoints', times(2:end))];
				dphi(:, jj) = dphi_temp;

				low_dphi = [false; abs(dphi_temp) < thresh; false];
				starts = ~low_dphi(1:end-2) & low_dphi(2:end-1);
				stops = low_dphi(2:end-1) & ~low_dphi(3:end); 
				if sum(starts)==0
					interval = nan; 
				else
					interval = max(times(stops) - times(starts));
				end
				perst(jj) = interval;
				dphi_mean(jj) = nanmean(abs(dphi_temp));

			end
			out.persistence = perst;
			out.dphi_mean = dphi_mean;
			out.dphi = dphi;
			out.phi_std = phi_std;
			out.phi_uw = phi_uw;
			
        end
        
		function out = summary_comparison(F, varargin)
			% out = summary_comparison(F, ::fig=gcf::, ::'save'::)
			% Generates stats comparing pairs of metrics
			SAVE = false;
			fig = gcf;
			for arg = varargin
				if isa(arg{:}, 'matlab.ui.Figure'), fig = arg{:};
				elseif strcmpi(arg{:}, 'save'), SAVE = true; 
				end
			end
			
			metrics = fieldnames(F);
			metrics(strcmpi(metrics, 'name')) = [];
			nF = numel(F.Name);

			metricpairs = nchoosek(metrics, 2);
			nM = size(metricpairs, 1);

			nrows = nM * nF;
			filename = cell(nrows, 1);
			whichpair = zeros(nrows, 1, 'uint16');
			dZ = cell(nrows, 1);
			[m1, R, theta, kappa, conf, sigma, N] = ...
				deal(nan(nrows, 1));

			idx = 0;
			for f = 1:nF % for each file


				name = F.Name{f};
			% 	data = load(name);

				for m = metricpairs'  % and each pair of metrics
					idx = idx + 1;

					whichpair(idx) = mod(idx-1, nM) + 1;
					filename{idx} = name;
					data1 = F.(m{1})(f);
					data2 = F.(m{2})(f);

					dd = data1.diff(data2);
					dZ{idx} = dd;

					N(idx) = numfinite_(dZ{idx});
					if N(idx) == 0, continue, end

			% 		m1(idx) = mean(dZ{idx}, 'omitnan');
			% 		R(idx) = abs(m1(idx));
			% 		kappa(idx) = circ_kappa(dd(finite));
					sigma(idx) = circ_std(dd, [], [], 'omitnan');
					[theta(idx), ul] = circ_mean(dd, [], [], 'omitnan');
					conf(idx) = ul - theta(idx);
				end

			end
			stats = table(filename, whichpair, R, theta, kappa, conf, sigma, N, ...
				m1, dZ);
			if SAVE, save('direction_stats', 'stats', 'metricpairs'); end
			
			out.stats = stats;
			out.metricpairs = metricpairs;
			out.fig = summary_stats(stats, metricpairs, fig);
		end

	end
	
end

%% Local functions
function n = numfinite_(X)
% n = sum(isfinite(X), dim='all')
n = sum(isfinite(X), 'all');
end

function [fits, metrics] = foldercontents2struct_(path, metrics)
    if nargin < 2 || isempty(metrics)
        metrics = {dir([path filesep '*.mat']).name};
        metrics = cellfun(@(ss) strrep(ss, '.mat', ''), ...
                metrics, 'uni', 0);
    end
    for ff = metrics
        fit = load([path filesep ff{:}]);
        pn = fieldnames(fit);
        fits.(ff{:}) = fit.(pn{:});
    end
end

function metrics = findsharedmetrics_(sz)
    N = numel(sz.patient);
    mtc = arrayfun(@(ii) ...
        dir(sprintf('WaveFits/%s_Seizure%d/*mat', ...
            sz.patient{ii}, sz.seizure(ii))...
        ), 1:N, 'uni', 0);
    temp = cellfun(@(xx) {xx.name}, mtc, 'uni', 0);
    [cc, M] = histcounts(categorical(cat(2, temp{:})));
    metrics = strrep(M(cc == N), '.mat', '');
        
end