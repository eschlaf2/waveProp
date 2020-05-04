classdef WaveProp
	
	properties 
		Name
		Vx = nan
		Vy = nan
		Magnitude = nan
		Direction = nan
		p = nan
		Beta = nan(1, 6)
		Curvature = nan
		Data
		t0
		sig = 0.05
		NClust
		ClustSize
		Quadratic = false
% 		HalfWin
% 		FBand
	end
	
	properties (Transient = true, Hidden = true)
% 		t_inds
% 		res
		
		Time
		Position
		scale_quiver = 1
		complexZ
		
	end
	
	properties (Hidden = true)
		WPParamNames = "Quadratic"
	end
	
	properties (Dependent = true, Hidden = true)
		time
		Z
		mask
		logp
		NormedMagnitude
	end
	
	methods 
		function ax = plot(obj, varargin)
			% ax = plot(type='2D', t0, ax); % if object has only one time point, t0=obj.t0
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
			% ax = plot(t0, ax); % if object has only one time point, t0=obj.t0
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
			if numel(obj.t0) == 1
				F = obj;
				units = '';
			else				
				[~, ind] = min((obj.t0 - tt).^2);
				F = obj.sub(ind);
				F.Magnitude = obj.NormedMagnitude(ind);
				units = 'std';
			end
			
			N = size(F.Data);
			[p1, p2] = ind2sub(N, find(isfinite(F.Data)));
			scatter3(ax, p1, p2, F.Data(isfinite(F.Data)), [], F.Data(isfinite(F.Data)), 'filled'); 
			hold(ax, 'on');
			[xx, yy] = ndgrid(1:N(1), 1:N(2));
			zfit = F.Beta(1) + F.Beta(2) * xx + F.Beta(3) * yy;
			im = surf(ax, xx, yy, min(F.Data(:)) * ones(N), F.Data, 'linestyle', 'none'); 
			im.Tag = 'Data';
			surfax = surf(ax, xx, yy, zfit, 'linestyle', 'none', 'facealpha', .5);
			surfax.Tag = 'Fit';
			title(ax, sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', F.t0, F.p, F.Magnitude, units)); 
			hold(ax, 'off');
			ax.Tag = checkname(['figs' filesep obj.Name '_' class(obj) '_' num2str(F.t0) '_plot3D']);
			
		end
		function mgn = get.NormedMagnitude(obj)
			mgn = obj.Magnitude;
			mgn(isoutlier(mgn)) = nan;
			mgn = (mgn - nanmean(mgn)) / nanstd(mgn);
		end
		function t0 = get.t0(obj)
			t0 = obj.t0(:);
		end
		function obj = set.Beta(obj, val)
			assert(isempty(val) || any(size(val) == 6));
			if isempty(val) || size(val, 2) == 6
				obj.Beta = val;
			else
				obj.Beta = val';
			end
		end
		
		function ax = plot2D(s, varargin)
			% ax = plot2D(tt, ax)
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
			if ndims(s.Data) == 3
				[~, which_t] = min(abs(tt - s.t0));
				data = s.Data(:, :, which_t);
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
				'MaxHeadSize', 1)
			hold(ax, 'off')
			ax.Tag = checkname(['figs' filesep obj.Name '_' class(obj) '_' num2str(F.t0) '_plot2D']);
		end
		
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
		
		function sref = subsref_nightmare(obj, s)
			
			switch s(1).type
				case '.'
% 					mtd = methods(obj, '-full');
% 					staticM = contains(mtd, 'Static') & contains(mtd, s(1).subs);
% 					if any(staticM)
% 						s(2).subs = [obj, s(2).subs];
% 					end
					sref = builtin('subsref', obj, s);
				case '()'
					if numel(obj) > 1
						sref = builtin('subsref', obj, s);
						return
					end
					sref = obj;
					N = numel(obj.t0);
					for f = string(fieldnames(obj)')
						temp = obj.(f);
						sz = size(temp);
						if strcmpi(f, 'data'), dim = 3; else, dim = 1; end
						
						if ~any(sz == N)
							sref.(f) = temp;
						elseif isvector(temp)
							sref.(f) = builtin('subsref', temp, s);
						else
							temp = reshape(shiftdim(temp, dim - 1), N, []);
							sz = circshift(sz, 1 - dim);
							data = temp(s.subs{1}, :);
							if numel(sz) < 3, sz = [sz 1]; end %#ok<AGROW>
							sref.(f) = shiftdim(reshape(data, sz(2:end)), numel(sz) - (dim - 1));
						end
					
					end
				case '{}'
					error('WaveProp:subsref',...
						'Not a supported subscripted reference')
				end
		end
% 		function pos = get.Position(obj)
% 			if isempty(obj.Position)
% 		end
		function sN = refit_data(s)
			W = warning;
			warning('off', 'stats:statrobustfit:IterationLimit');
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
			N = numel(s.t0);
			s.Beta = nan(N, 6);
			s.Curvature = nan(N, 1);
			
			dims = size(s.Data);
			[xx, yy] = ndgrid(1:dims(1), 1:dims(2));
			position = [xx(:) yy(:)];
			for ii = 1:N
				data = s.Data(:, :, ii);
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
		end
		function [H, centers, T] = hist(obj, window, mean_center, show_directions, h)
			% [H, edges, T] = hist(window=1, mean_center=true, show_directions=true, h=gcf)
			
			if nargin < 5, h = gcf; end
			if nargin < 4 || isempty(show_directions), show_directions = true; end
			if nargin < 3 || isempty(mean_center), mean_center = true; end
			if nargin < 2 || isempty(window), window = 1; end
			t = obj.t0;
			if mean_center
				obj.Direction = angle(exp(1j*(obj.Direction - obj.mean)));
			end
% 			dir = exp(obj.Vx + 1j*obj.Vy);
			T = t(1):.01:t(end)+.01;
			edges = linspace(-pi, pi, 65);
			H = nan(length(edges)-1, length(T));
			for ii = 1:length(T)
				t_inds = abs(t - T(ii)) <= window/2;
				dir = obj.Direction(t_inds);
				H(:, ii) = histcounts(dir, edges) ./ window;
			end
			set(0, 'currentfigure', h);
			ax = subplot(1, 9, 1:6);
			s = subplot(1, 9, 7:8);
			win = gausswin(round(pi/4 / diff(edges(1:2)))) .* gausswin(round(window ./ diff(T(1:2))))';
			win = win / sum(win, 'all');
			H = [H; H; H];
			centers = unwrap(repmat(edges(1:end-1)' + diff(edges(1:2)/2), 3, 1)) - 2*pi;
			H = conv2(H, win, 'same');
			angle_mask = abs(centers) <= pi;
			centers = centers(angle_mask);
			H = H(angle_mask, :);
			contourf(ax, T, centers, H, quantile(H(:), .1:.05:.9), 'linestyle', 'none'); %conv2(H, ones(4, window/2*1e3)/(2*window*1e3), 'same'));
			if show_directions
				NP = ax.NextPlot;
				ax.NextPlot = 'add';
				plot(ax, obj.t0, obj.Direction, '.', 'markersize', 10);
				ax.NextPlot = NP;
			end
			C = colorbar(ax);
			temp = subplot(1, 9, 9); 
			temp.Visible = 'off';
			C.Position(1:2) = temp.Position(1:2);
% 			axis(ax, 'xy');
			plot(s, smoothdata(max(H, [], 2)), centers, 'linewidth', 2);
			axis(s, 'tight');
			xticks(s, []); yticks(s, []);
			xlabel(ax, 'Time (s)')
			ylabel(ax, 'Direction')
% 			SR = 1/min(diff(t));
% 			W = ones(floor(SR * window), 1);
% 			T = min(t):1/SR:max(t)+1/SR;
% 			D = zeros(size(T));
% 			t_inds = interp1(T, 1:length(T), t, 'nearest');
% 			D(t_inds) = dir;
		end
		function logp = get.logp(s)
			logp = -log(s.p)/log(5);
		end
		function Z = get.complexZ(s)
			Z = complex(s.Vx, s.Vy);
		end
		function mask = get.mask(s)
			mask = s.p > s.sig | isnan(s.Magnitude);
		end
		function time = get.time(s)
			time = s.t0;
		end
		function Z = get.Z(s)
			Z = s.Direction;
		end
		function Vx = get.Vx(s)
			Vx = s.Vx;
			Vx(s.mask) = nan;
		end
		function Vy = get.Vy(s)
			Vy = s.Vy;
			Vy(s.mask) = nan;
		end
		function Curvature = get.Curvature(s)
			Curvature = s.Curvature;
			Curvature(s.mask) = nan;
		end
		function D = get.Direction(s)
			D = s.Direction;
			D(s.mask) = nan;
		end
		function M = get.Magnitude(s)
			M = s.Magnitude / 400;  % Fits return electrodes/second
			M(M > 1e6) = nan;
			M(isoutlier(M)) = nan;
% 			M(M./nanstd(M) > 4) = nan;
		end
		
		function [mn, ci] = mean(s, data)
			if nargin < 2, data = s.Direction; end
			[mn, ul] = circ_mean(data, [], [], 'omitnan');
			ci = ul - mn;
		end
		
		function F = sub(obj, ind)
			N = min(numel(obj.t0), numel(obj.Vx));
			eval(sprintf('F = %s;', class(obj)));
			for f = string(fieldnames(obj)')
				temp = obj.(f);
				sz = size(temp);

				if all(sz < N) || ischar(obj.(f))
					F.(f) = temp;
				elseif isvector(temp)
					F.(f) = temp(ind);
				elseif ismatrix(temp)
					F.(f) = temp(ind, :);
				else
					assert(size(temp, 3) >= N);
					dim = 3;
					temp = shiftdim(temp, dim-1);
					sz = size(temp);
					temp = reshape(temp, N, []);
					data = reshape(temp(ind, :), sz(2:3));
					F.(f) = data;
				end

			end
		end
		
		function F = split(obj)
			N = min(numel(obj.t0), numel(obj.Vx));
			F(N) = obj.sub(N);
			for ii = 1:N
				F(ii) = obj.sub(ii);
			end
		end
		
		function d = diff(s, other)
			t_inds = interp1(other.t0, 1:length(other.t0), s.t0, 'nearest', 'extrap');
			d = angle(exp(1j*(s.Direction - other.Direction(t_inds))));
		end
		
		function [d, xi, bw] = ksdensity(obj, ax, rotateby, xres, bw, varargin)
			% [d, xi, bw] = ksdensity(ax=gca, rotateby=0, xres=pi/500, bw=.15*pi/3, ::plot directives::)
			if nargin < 5 || isempty(bw), bw = .15*pi/3; end
			if nargin < 4 || isempty(xres), xres = 500; end
			if nargin < 3 || isempty(rotateby), rotateby = 0; end
			if nargin < 2 || isempty(ax), ax = gca; end
			temp = angle(exp(1j * (obj.Direction - rotateby)));
			temp = [temp(:) - 2*pi; temp(:); temp(:) + 2*pi];
		% 	ksdensity([repmat(ct.(f{:}), 3, 1) temp], [x1 x2], 'bandwidth', [5 .15*pi/3]);
			gridx2 = linspace(-3*pi, 3*pi, 6*xres + 1);
			[d, xi, bw] = ksdensity(temp, gridx2, 'bandwidth', bw);
			plot(ax, xi, d/max(d), varargin{:});
			xlim(ax, [-pi pi])
		end
	
		
		function S = compile(obj, t_fun)
			% Compile an array of WaveProp objects into one
			% S = obj.compile(t_fun=@(t) t0 < Inf)
			if nargin < 2, t_fun =@(t0) t0 < Inf; end
			S = obj(1);
			for ii = 1:numel(obj)
				t_mask = t_fun(obj(ii).t0);
				obj(ii).p(t_mask) = Inf;
			end
% 			S.p = cat(1, obj.p);
% 			S.Magnitude = cat(1, obj.Magnitude);
			for ff = string(fieldnames(obj)')
				
% 				if ismember(lower(ff), {'name'}), continue, end
				try
					S.(ff) = cat(1, obj.(ff));
				catch ME
					S.(ff) = {obj.(ff)};
					if ~strcmpi(ME.identifier, 'MATLAB:catenate:dimensionMismatch')
						rethrow(ME)
					end
				end
			end
% 			t_mask = t_fun(S.t0);
% 			S.p(t_mask) = Inf;
			t0d = [0; diff(S.t0)];
			t0d(t0d <= 0) = 1;
			S.t0 = cumsum(t0d);
			S.sig = min(S.sig);
			
		end
		
		function obj = parse_inputs(obj, varargin)
			for ii = 1:2:numel(varargin)
				ff = validatestring(varargin{ii}, [obj.ParamNames(:); obj.WPParamNames(:)]);
				obj.(ff) = varargin{ii+1};
			end
		end
	end
	
%% Static methods
	
	methods (Static)
		function [ax, tt, directive] = parse_plot_inputs(varargin)
			ax = gca; tt = 0;
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
		end
		
		function [mat, tt] = WP2mat(S, field)
			N = numel(fieldnames(S));
			tt = [];
			for ff = string(fieldnames(S)')
				tt = [tt; S.(ff).t0(:)]; %#ok<AGROW>
			end
			tq = unique(tt, 'stable');
			X = {};
			for ff = string(fieldnames(S)')
				[tt, ii] = unique(S.(ff).t0);
				iq = interp1(tt, ii, tq, 'nearest', 'extrap');
				X = [X {S.(ff).(field)(iq)}]; %#ok<AGROW>
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

			MIN_RATIO_FINITE = 0.25;                                                 % we require 50% of electrodes to have a defined delay

			P0 = nan;
			b = nan;
			
			finite = numfinite(data);
			
			
			if finite > 3 && finite/numel(data) >= MIN_RATIO_FINITE  % check enough delay data is not NaN.
				if any(ismember({'quadratic', 'quad'}, lower(varargin)))
					X = [position prod(position, 2) position.^2];  % predictors
				else 
					X = position;
				end
				[b,stats] = robustfit(X, data, 'fair');                     % fit the delay vs two-dimensional positions
% 				H = [0 1 0; 0 0 1];  % These are default values
% 				c = [0 ; 0];
				P0 = linhyptest(b, stats.covb, [], [], stats.dfe);                    % perform F test that last two coefs are both 0.
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
			S.Magnitude = cat(1, s.Magnitude);
			for f = string(fieldnames(s)')
				if ismember(lower(f), {'quadratic', 'sig', 'halfwin', 'fband', 'name', 'minfreq', 'p', 'magnitude'}), continue, end
				if isempty(S.(f)), continue,
				elseif ismember(lower(f), {'data', 'freq'})
					S.(f) = cat(3, s.(f));
				else
					S.(f) = cat(1, s.(f));					
				end
			end
		end
		function obj = struct2obj(fstruct)
			% Usage: obj = struct2obj(res)
			
			obj = WaveProp;
			assert(nargin == 1);
			obj.Name = fstruct.Name;
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
			
			% Cleaning
			if args.original
				fname =@(p, s) sprintf('%s_Seizure%s_Neuroport_10_10_wave_prop.mat', p, s);
			else
				fname =@(p, s) strrep(sprintf('%s*Seizure*%s*fits.mat', p, s), '***', '*');
			end
			if isempty(args.files)
				if ischar(args.pat) && ischar(args.seizure)
					if strcmpi([args.pat args.seizure], '**')
						fid = fopen('seizures2.txt');
						A = textscan(fid, '%s %s'); 
						pats = A{1}; seizures = A{2};
						for ii = numel(pats):-1:1
							files(ii) = dir(fname(pats{ii}, seizures{ii}));
						end
						fclose(fid);
						args.files = files;
					else
						args.files = dir(fname(args.pat, args.seizure));
					end
				else

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
		function fits = load(varargin)
			s = [];
			if nargin > 0 && isa(varargin{1}, 'MEA')
				s = varargin{1};
				varargin(1) = [];
			end
			args = WaveProp.parse_load(varargin{:});
			if isa(s, 'MEA')
				if args.original
					fname = sprintf( ...
						'%s_Seizure%s_Neuroport_10_10_wave_prop.mat', ...
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
					if ii == nF && isempty(metrics)
						res = load([files(ii).folder filesep files(ii).name]);
						metrics = fieldnames(res); 
					else
						res = load([files(ii).folder filesep files(ii).name], metrics{:});
						not_member = ~ismember(metrics, fieldnames(res));
						if any(not_member)
							fits = rmfield(fits, metrics(not_member));
							metrics(not_member) = [];
						end
					end
					if args.original
						fits(ii) = WaveProp.file2fit(res);
					else
						res.Name = strrep(files(ii).name, '.mat', '');
						fits(ii) = res;
					end
				end	
				for ff = string(fieldnames(fits)')
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
				dir_stats{pp} = D;
				
				for mm = 1:numel(metrics)
					if numfinite(magnitude(:, mm)) < 10, continue; end
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
		function out = compute_dir_stats(phi, t0, thresh)
			if nargin < 3, thresh = pi/64; end
			[dphi, phi_uw, phi_std] = deal(nan(size(phi)));
			for jj = size(phi, 2):-1:1
				phi_temp = phi(:, jj);  % get Z for given patient and metric

				mask = isfinite(phi_temp);  % isolate non-nan
				phi_uw_temp = unwrap(phi_temp(mask));  % unwrap (isolate and then unwrap is questionable ...)
				phi_temp(mask) = phi_uw_temp;
				phi_uw(:, jj) = phi_temp;
				for tt = 1:numel(t0)
					mask_tt = (t0 >= t0(tt) - .5) & (t0 <= t0(tt) + .5);
					phi_std(tt, jj) = nanstd(phi_temp(mask_tt));
				end
		% 		phi_temp_sm = smoothdata(phi_temp, 1, 'gaussian', 1, 'SamplePoints', t_temp);  % smooth with a 1 second (gaussian) kernel
		% 		dphi{ii}(:, jj) = [nan; smoothdata(diff(phi_temp_sm))];
				dphi_temp = [nan; smoothdata(diff(phi_temp), 1, 'movmean', 1, 'samplepoints', t0(2:end))];
				dphi(:, jj) = dphi_temp;

				low_dphi = [false; abs(dphi_temp) < thresh; false];
				starts = ~low_dphi(1:end-2) & low_dphi(2:end-1);
				stops = low_dphi(2:end-1) & ~low_dphi(3:end); 
				if sum(starts)==0
					interval = nan; 
				else
					interval = max(t0(stops) - t0(starts));
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
			% out = summary_stats(F, ::fig=gcf::, ::'save'::)
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

					N(idx) = numfinite(dZ{idx});
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