classdef WaveProp < handle & matlab.mixin.Copyable
	% I want this to take a set of TOAs and show me things about them... So
	% then, I think I want my MaxDescent and Delays to only hold the data
	% and method-specific parameters (this is probably not far off of what
	% they do)
    
    
    %%% Commenting these out (5/10/21); may need to make other adjustments
    %%% as a result, but I'm pretty sure these are things that just reflect
    %%% inefficiencies of learning
%     properties (Hidden = true)
% 		WPParamNames = "Quadratic"
% 	end

	properties 
		
        Patient string
        Seizure string
        Position
        UseLargestCluster = true
        
        GridSize = [10 10]  % Size of the MEA (# of electrodes along each dimension; always [10 10] unless it's a sim)
		Length = 4  % length of the MEA (always 4 mm, unless it's a sim)
        MMpElectrode = .4  % mm per electrode (again, always 4 unless it's a sim)
        
        TOA  % this is the important one to save (especially for delays method since it takes a long time to compute
        TOAcube  % converts TOA into a TxXxY array
        t0

    end
    

	properties (Transient = true)

        ResultOfFit
        NClust  % I don't even know what this does anymore (5/9/21)
		ClustSize
        
        sig = 0.05
        RotateBy = 0
        MinFinite = 30
        Quadratic = false
		Original = false
		Magnitude
        AltMagnitude
		Direction
		
		scale_quiver = 1
		complexZ
		Inds = []
		Early = [-Inf 25]
		Late = [30 Inf]
		MinDetections = 100
        
        % Move these to dependent after the next rerun (of everything...)
        Vx 
		Vy
		p
		Beta
		Curvature
	end
	
	
	
	properties (Dependent = true, Hidden = false)

        Name string
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
        Data
        locs
    end
    
    methods  % imports
        [H, centers, T] = hist(obj, window, mean_center, show_directions, h)
        ax = particle_smoother(W, ax, NPARTS, SPREAD_FACTOR, SD_WIN, ANG_RES, SMOOTHING_WIN, SHOW_DIRECTIONS, SHOW_EVOLUTION)
    end
	
	methods  % Getters and Setters
        
        function set.Name(obj, value)
            value = string(value);
            ps = strsplit(value, {'_Seizure', ' ', '_'});  % Assume PP_SeizureXX assignment
            obj.Patient = ps(1);
            obj.Seizure = ps(2);
        end
        
        function name = get.Name(obj)
            name = sprintf("%s_Seizure%s", obj.Patient, obj.Seizure);
        end
        
        function [ci_dir, ci_sp] = direction_ci(M)
            % Estimate a confidence interval around the estimated direction
            % by bootstrapping (as in Martinet, 2017)
            
            res = M.ResultOfFit;
            N = size(res, 1);
            
            [ci_dir, ci_sp] = deal(nan(N, 2));
            
            for ii = 1:N
                if isnan(res.stats(ii).covb(1)), continue; end
                samples = mvnrnd(res.beta(ii, 2:3), res.stats(ii).covb(2:3,2:3), 1000);       % generate samples from model.
                boot_dir = atan2(samples(:,2), samples(:,1));                   % bootstrap direction
                boot_sp = 1./sqrt(samples(:,1).^2 + samples(:,2).^2);           % bootstrap velocity.

                mn = circ_mean(boot_dir);
                ci_dir(ii, :) = ...
                    quantile( ...
                        angle(exp(1j * (boot_dir - mn))), ...
                    [0.025, 0.975]) ...
                    + mn;                    % CI for direction
                ci_sp(ii, :) = quantile(boot_sp, [0.025, 0.975]);  
            end
            ci_dir = angle(exp(1j * (ci_dir - M.RotateBy)));
        end
        
        function out = get.ResultOfFit(obj)
            % If fits haven't been run, run them now (fit a plane to each
            % time)
            
            if isempty(obj.ResultOfFit)
                inds_ = obj.Inds; % Temporarily set inds to default to only generate fits for actual data
                obj.Inds = [];
                
                
                % Get the observed TOA
                if obj.UseLargestCluster
                    toa = obj.use_largest_cluster(obj.TOAcube);
                    toa = toa(:, obj.locs);
                else
                    toa = obj.TOA;
                end
                Nt = size(obj.TOA, 1);
                
                % Initialize matrices to store the results and fit a plane
                % at each time point
                V = nan(Nt, 2); p_ = nan(Nt, 1); beta = nan(Nt, 3);
                stats_temp = cell(Nt, 1);
                for ii = 1:Nt
                    [V(ii, :), p_(ii), beta(ii, :), stats_temp{ii}] = ...
                        obj.fit_plane(toa(ii, :), obj.Position);
                end
                
                % Clean up stats object
                mask_ = cellfun(@isstruct, stats_temp);
                tpl = stats_temp{find(mask_, 1)};
                for ff = string(fields(tpl))', tpl.(ff)(:) = nan; end
                stats_temp(~mask_) = {tpl};
                stats = cat(1, stats_temp{:});
                
                % Save the result to a table
                obj.ResultOfFit = table(V, p_, beta, stats);
                
                % Reset inds to whatever they were
                obj.Inds = inds_;
            end
            out = obj.ResultOfFit;
        end
        
        function toa = get.TOA(obj)
            if ndims(obj.TOA) == 3
                obj.TOA = obj.TOA(:, obj.locs);
            end
            toa = reshape(obj.TOA, numel(obj.t0), size(obj.Position, 1));
            toa = toa(obj.Inds, :);
        end
        
        function res = get.TOAcube(obj)
            % Convert TOA to 3D (and return .Inds)
            
            m = numel(obj.time); 
            toa_size = [m, obj.GridSize];  % result should be TIMExXxY
            
            res = nan(toa_size);  % initialize matrix to hold the result

            % Get the indices where the TOA data belongs and generate the
            % resulting TOA
            [tt, xx] = ndgrid(1:m, obj.Position(:, 1));
            [~, yy] = ndgrid(1:m, obj.Position(:, 2));

            locs_ = sub2ind(toa_size, tt(:), xx(:), yy(:));
            res(locs_) = obj.TOA;

            

        end
        
        function data = get.Data(self)
            % Wrapper for self.TOAcube
            % renamed Data to TOA to be more clear about what's in here
            
            data = self.TOAcube;
        end
        
        function set.Data(self, value)
            % Data is a wrapper for TOA (more clear naming)
            self.TOA = value;
        end
        
        function fname = get.FName(self)
            % Returns the name of the file where the TOA are stored
            fname = sprintf("WaveFits/%s_Seizure%s/", self.Patient, self.Seizure);
            warning('You probably want to call .Name here rather than .FName');
%             fname = sprintf('%s_Seizure%s_fits.mat', self.Patient, self.Seizure);
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
        
        function set.time(self, value)
            self.t0 = value;
        end
        
		function beta = get.Beta(obj)
			beta = obj.ResultOfFit.beta(obj.Inds, :);
        end
        
        function vx = get.Vx(obj)
            vx = obj.ResultOfFit.V(:, 1);
        end
        
        function vy = get.Vy(obj)
            vy = obj.ResultOfFit.V(:, 2);
        end
        
		function mgn = get.NormedMagnitude(obj)
			mgn = filloutliers(obj.Magnitude, nan, 'median', 'thresholdfactor', 10);
			mgn = normalize(mgn);
			mgn = mgn(obj.Inds);
        end
        
		function logp = get.logp(s)
            % returns XX from expression p=5e-XX
			logp = -log10(s.p/5);
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
			M = filloutliers(sqrt(s.Vx.^2 + s.Vy.^2), nan, 'median', ...
                'thresholdfactor', 100);
            
            % require full window (set to .95 bc a few of the windows are short a ms)
            if ismember('T', properties(s))  
                M(s.T < .95*2*s.HalfWin) = nan;
            end
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
        
		function D = get.Direction(M)
            
            if isempty(M.Direction)  % this should usually be empty so it pulls from (Vx, Vy) but sometimes I want WaveProp functions on a specific set of directions
                D = atan2(M.Vy, M.Vx);

                inds_ = M.Inds;
                M.Inds = [];
                D(M.mask) = nan;

                % Remove points where directions show large shifts within 50
                % ms. (this shouldn't do anything to M and D10 methods, but some
                % experimental methods [e.g. M10] benefit from this)
                dir_sm = movmean(exp(1j*D), .05, 'samplepoints', M.t0);  
                dir_sm(abs(dir_sm) < cos(pi/8)) = nan; % this ensures consecutive angles differ by less than 45°
                d2 = movmean(dir_sm, .04, 'omitnan', 'samplepoints', M.time);  % waves traveling at 100mm/s should be on the MEA ~40 ms
                d2(abs(d2) < .9) = nan;  % require strong similarity in directions


                D = angle(d2) - M.RotateBy;
                D = angle(exp(1j * D));  % Keep in range [-pi pi]
                M.Inds = inds_; 
                D = D(inds_);
            else
                D = M.Direction;
            end
        end
        
		function M = get.Magnitude(self)
            % in mm/s
			M = sqrt(self.Vx.^2 + self.Vy.^2) * self.MMpElectrode;  % Fits return electrodes/second. Convert to mm/s
			M = M(self.Inds);
            M(self.mask) = nan;
        end
        
        function M = get.AltMagnitude(self)
            t_ext = range(self.Data, [2 3]);  % time taken for wave to cross MEA
            M = self.Length * sqrt(2) ./ t_ext;  % speed in mm/second (5.67 mm = length of MEA diagonal)
            M = M(self.Inds);
            M(self.mask) = nan;
        end
        
		function psig = get.p(s)
			psig = s.ResultOfFit.p_(s.Inds);
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
	
	methods  % Most of the things
		
        %% Helpful others
        
        function ll = get.locs(M)
            % convert x-/y-position to indices
            ll = sub2ind(M.GridSize, M.Position(:, 1), M.Position(:, 2));
        end
        function t_inds = time2inds(M, tt)
            [~, t_inds] = min(abs(M.time(:) - tt(:)'));
        end
        
        function [dat, tt, pos] = preproc_discharges(M)
            
            % Reshape data into a 2D matrix; 
            BASELINE_CUTOFF = 2;  % Use 4 if you switch to robust normalization
            dat = reshape(M.Data, length(M.Data), []);
            tt = M.time;
            
            % Exclude channels that are rarely active
            locs_ = find(normalize(sum(isfinite(dat))) > -1);
            dat = dat(:, locs_);
            
            % If using very short windows (i.e. likely to have only a few
            % electrodes active at a time), "extend" the windows using a
            % moving average
            if M.HalfWin < .03
                dat = movmean(dat + tt, .1, 'omitnan', 'SamplePoints', tt);
                dat = filloutliers(dat, nan, 2);  % remove outliers introduced by smoothing
            end
            
            % Remove time points with fewer than MinFinite electrodes
            % active
            mask_ = sum(isfinite(dat), 2) < M.MinFinite;
            dat(mask_, :) = [];
            tt(mask_) = [];
            
            % Require that discharges be at least 50 ms apart (again, this
            % shouldn't really be doing anything in the M or D10 methods
            % since those are already controlling for this... mostly for
            % experimental short window methods)
            
            % Reshape the data so there are zeros between discharge times
            ttI = tt(1):.001:tt(end);
            [~, tt_inds] = min(abs(ttI - tt), [], 2);
            countsI = zeros(size(ttI));
            countsI(tt_inds) = sum(isfinite(dat), 2);
            
            [~, locs_t] = findpeaks(countsI, ttI, ...
                'MinPeakDistance', .05, 'minpeakheight', M.MinFinite);
            [~, locs_i] = min(abs(locs_t - tt));
            
            % Z-score the data at each time point
            dat = normalize(dat(locs_i, :), 2, 'zscore', 'std');
            tt = tt(locs_i);
            nanmask = all(isnan(dat), 2);  % robust normalizing will remove all points sometimes (shouldn't need it after prior treatment anyway)
            dat(nanmask, :) = [];
            tt(nanmask) = [];
            
            % Remove forms that look nothing like anything else (i.e.
            % discharges that have high distance at the .01 quantile)
            ddd = pdist(fillmissing(dat, 'constant', 0), 'squ');
            ddd = squareform(ddd/size(dat, 2));
            
            level = quantile(ddd, .01);
            
            locs_i = level < BASELINE_CUTOFF;
            if sum(~locs_i) > 2
                % Again, this shouldn't be doing much. It's here for
                % the sake of the clustering more than anything, but noise
                % should be well removed by now, so send a warning if this
                % removes more than two discharges. 
                warning('More than 2 discharges (%d) removed as outliers.', sum(~locs_i)); 
            end
            
            dat = dat(locs_i, :);
            tt = tt(locs_i);
            
            % Compute position for easy shaping on return
            [xx, yy] = ind2sub([10, 10], locs_);
            pos = [xx(:), yy(:)];
            
        end
        
        function [dat, tt, pos, IDX, C, D] = cluster(M, h, varargin)
            if nargin < 2 || isempty(h)
                h = figure; fullwidth(true); 
            end
            T = tiledlayout(h, 3, 5); 
            DIST = 'corr';
            
            k = 1;
            
            [dat, tt, pos] = M.preproc_discharges;
                        
            X = fillmissing(dat, 'constant', 0);
            % eva = evalclusters(single(X), @kmedoids, 'gap', 'klist', 5:30);

            N = size(dat, 1);
            klist = 2:10;

            W = warning; warning('off');
            switch DIST
                case 'sqeuc'
                    dist_fun = @(x, k) kmeans(x, k, 'distance', 'sqeuc', ...
                        'replicates', 500, 'emptyaction', 'drop');
                    drop_fun = @(D) min(D, [], 2) > size(dat, 2);
                case 'corr'
                    dist_fun = @(x, k) kmeans(x, k, 'distance', 'corr', ...
                        'replicates', 500, 'emptyaction', 'drop');
                    drop_fun = @(D) min(D, [], 2) > .5;
            end
            % eva = evalclusters(X, dist_fun, 'davies', 'klist', klist);


            [IDX, ~, ~, D] = dist_fun(X, floor(N/10));  % Get clusters

            % Drop small clusters and bad matches
            cN = arrayfun(@(x) sum(x == IDX), IDX);
            drop = drop_fun(D) | cN(:) < 3;
            IDX(drop) = -1;

            % Cluster again without dropped discharges
            eva2 = evalclusters(X(IDX > 0, :), dist_fun, 'davies', 'klist', klist);
            [IDX2, C2, ~, D2] = dist_fun(X(IDX > 0, :), eva2.OptimalK);
            chg_id = IDX > 0;
            IDX(chg_id) = IDX2;
            C = C2;
            D = nan(N, max(IDX2));
            D(chg_id, :) = D2;

            % Drop again
            cN = arrayfun(@(x) sum(x == IDX), IDX);
            drop = drop_fun(D) | cN(:) < 10;
            IDX(drop) = -1;

            % Re-number
            [G, bigC] = findgroups(IDX);
            bigC(bigC == -1) = [];
            C = C(bigC, :);
            D = D(:, bigC);
            IDX = G - 1;


            nexttile(T, k); k = k + 1;
            % imagesc(squareform(pdist(C, 'squ')/size(dat, 2)), [0 2]); colorbar; axis square
            imagesc(squareform(pdist(C, 'corr')), [0 2]); colorbar; axis square
            cmap = make_diverging_colormap([1 0 0; 1 1 1; 0 0 0]);
            colormap(gca, cmap)
            title('pdist clust')


            nexttile(T, k, [1, 2]); k = k + 2;
            plot(tt, IDX, '.')
            ylim([-1 max(IDX) + 1])
            title(sprintf('Clust id v. time (N=%d, %0.2f)', N, sum(IDX > 0)/N))


%             [xx, yy] = ind2sub([10, 10], locs_);
            xx = pos(:, 1); yy = pos(:, 2);
            cmap = make_diverging_colormap([1 0 0; 0 0 0], .85*[1 1 1]); 
            k = max(k, 6);
            for ii = 1:size(C, 1)
                nexttile(T, k); k = k+1;
                cN = sum(IDX == ii);
                scatter(xx, yy, 100, C(ii, :), 'filled');
                xlim([0 11]); ylim([0 11]);
                title(sprintf('%d: %d (%0.2f)', ii, cN, cN/N))
                axis square
                colormap(gca, cmap)
            end

            warning(W)
            title(T, strrep(M.Name, '_', ' '));
            
        end
        function D2 = nan_euc_dist(M, XI, XJ)
            % Don't think this is used anywhere...
            n = size(XI,2);
            sqdx = (XI-XJ).^2;
            nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
            nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
            D2squared = nansum(sqdx,2).*n./nstar; % Correction for missing coordinates
            D2 = sqrt(D2squared);
        end
        
        function [dist_rho, tpl_id, tt, full] = distance(M, templates)
            % Converts correlation to distance (1 - rho);
            % [dist_rho, tpl_id, tt, full] = distance(M, templates)
            [rho, ~, tt, rho_sig, pval_sig] = M.correlation(templates);
            dist_rho = 1 - rho_sig;
            full = 1 - rho;
            tpl_id = pval_sig(:, 2);
            
        end
        
        function [rho, pval, tt, rho_sig, pval_sig] = correlation(M, templates)
            % [rho, pval, tt, rho_sig, pval_sig] = correlation(M, templates)
            % Outputs:
            %     rho: correlation for each discharge-template pair. 
            %     pval: p-value (non-zero rho) corresponding to rho
            %     tt: is the time of the discharge. 
            %     rho_sig: correlation of matching the lowest pvalue (most
            %         certainty of non-zeros correlation)
            %     pval_sig: (Nx2), lowest pvalue in first column and template
            %         number in second column
            
            assert(~isempty(templates), 'Input must be non-empty.')
            
            % z-score the templates and reshape to 2d (time x electrode#)
            templates = normalize(reshape(templates, size(templates, 1), 100)');

            % preprocess discharges (a little cleaning, but mostly this
            % just extracts the times where there are enough finite
            % electrodes and normalizes the data at each time point)
            [dat, tt, pos] = M.preproc_discharges;  
            locs_ = sub2ind([10 10], pos(:, 1), pos(:, 2));
            
            % initialize matrices to hold results
            rho = nan(size(dat, 1), size(templates, 2));
            pval = rho;
            
            % compute correlate (rho) and p-val (where H0 is "correlation
            % is 0") for each TW and template pair
            for ii = 1:size(templates, 2)
                for jj = 1:size(dat, 1)

                    % reshape the discharge to match the template shape and
                    % put template and discharge into a 100x2 array
                    temp = nan(100, 1);  
                    temp(locs_) = dat(jj, :); 
                    XY = [templates(:, ii), temp];
                    
                    % Remove electrodes with nan values
                    XY = XY(~any(isnan(XY), 2), :);
                    
                    % if fewer than 5 electrodes remain, don't compute corr
                    if size(XY, 1) < 5, continue; end
                    [rho_, p_] = corr(XY);
                    
                    % store the results
                    rho(jj, ii) = rho_(2);
                    pval(jj, ii) = p_(2);
                end
            end
            
            % Additionally return the correlation and p-val for the
            % template where each discharge had the most significant p-val.
            
            % Choose the template with the most significant pvalue
            [~, template_idx] = min(pval, [], 2);
            
            inds = sub2ind(size(rho), 1:size(rho, 1), template_idx');
            pval_sig = pval(inds(:));
            mask_ = isnan(pval_sig);
            rho_sig = rho(inds(:));
            rho_sig(mask_) = nan;
            template_idx(mask_) = nan;
            pval_sig = [pval_sig template_idx];
        end
        
        function [dir_out, sp_out, dir_sm, win, thresh] = discharge_directions(M, win, thresh)
            % out = discharge_directions(M, win=0.05, thresh=pi/4)
            % Returns the direction from each independent discharge.
            % Inputs: 
            %   win: min time between discharges to use (in seconds)
            %   thresh: threshold on directionality  
            if nargin < 2 || isempty(win), win = .05; end
            if nargin < 3 || isempty(thresh), thresh = pi/4; end
            
            D = M.Direction;
            tt = M.time;
            SP = M.AltMagnitude;
            SP(~isfinite(D)) = nan;
            
            dir = exp(1j * D);
            dir_sm = movmean(dir, win, 'omitnan', 'SamplePoints', tt);
%             dir_sm = movmean(dir, 2, 'omitnan');
            [~, locs_] = findpeaks( ...
                fillmissing(abs(dir_sm), 'constant', 0), tt, ...
                'minpeakheight', max(cos(thresh/2), .5), ...
                'minpeakdistance', win);
            [~, locs_i] = min(abs(tt(:) - locs_(:)'));
            dir_out = deal(nan(size(dir_sm)));
            dir_out(locs_i) = angle(dir_sm(locs_i));
            
            sp_out = movmean(SP, win, 'omitnan', 'SamplePoints', tt);
            sp_out(~isfinite(dir_out)) = nan;
            
        end
        
        function [mn_est, ci_est] = mov_mean_dir(M, win, times)
            if nargin < 3, times = M.time; end
            if nargin < 2 || isempty(win), win = 10; end
            
            WNG = warning;
            warning('off', 'circ_confmean:requirementsNotMet');
            
            dir = M.Direction;
            mask_ = isfinite(dir);
            
            mn_est = circ_mov_stat('mean', dir(mask_), times, win, ...
                'samplepoints', M.time(mask_));
            ci_est = circ_mov_stat('confmean', dir(mask_), times, win, ...
                'samplepoints', M.time(mask_));
            
            warning(WNG);
        end
        
        function dir_sm = smooth_angular(M, dir, win, thresh)
            % Returns the directions smoothed over <win> second windows
            % where directionality is at least <thresh> (i.e. length of
            % summed vectors is at least <thresh>)
            if nargin < 3 || isempty(win), win = .05; end
            if nargin < 4 || isempty(thresh), thresh = 0.7; end
            dir = exp(1j*dir);
            dir_sm = movmean(dir, win, 'omitnan', 'SamplePoints', M.time);
            dir_index = abs(dir_sm);
            dir_sm(isnan(dir) | dir_index < thresh) = nan;
            dir_sm = angle(dir_sm);
        end
        
        function [res, dir_sm, win, thresh] = ZZdischarge_directions(M, win, thresh)
            % out = discharge_directions(M, win=0.1, thresh=pi/4)
            % Returns the first direction from each independent discharge.
            % Inputs: 
            %   win: smoothing window to use (in seconds)
            %   thresh: threshold on directionality and allowable 
            %       difference between consecutive directions 
            
            if nargin < 2 || isempty(win), win = .05; end
            if nargin < 3 || isempty(thresh), thresh = pi/4; end
            
            D = M.Direction;
            tt = M.time;
            
            dir = exp(1j * D);
            locs_actual = find(isfinite(dir));
            dir_sm = movmean(dir, win, 'omitnan', 'SamplePoints', tt);
            dir_sm(abs(dir_sm) < max(cos(thresh/2), .5)) = nan;
            
%             dir_sm(isnan(movmean(dir, 2, 'omitnan'))) = nan;  % allow expansion by only 1
%             d2 = fillmissing(dir_sm, 'next');
%             d2(isnan(dir_sm)) = d2(isnan(dir_sm)) * -1;  % Flip by 180° so difference register
            
            
            diffs = diff(angle(dir_sm(:)));
            dir_diffs = angle(exp(1j*diffs));
            locsD = find(abs(dir_diffs) > thresh);
            N = movsum(isfinite(dir_sm), win, 'omitnan', 'samplepoints', tt);
            [~, locsN] = findpeaks(N);
            locs = unique([locsD(:)+1; locsN(:)]);
            ll = interp1(locs_actual, 1:length(locs_actual), locs, 'nearest', 'extrap');
            locs_t = locs_actual(ll);
            [G, locs_t] = findgroups(locs_t);  % in case two locs map to the same original
            dd = splitapply(@nanmean, dir_sm(locs), G);
            res = nan(size(dir_sm));
%             out(locs_t) = angle(dir_sm(locs));
            res(locs_t) = angle(dd);
            dir_sm = angle(dir_sm);
            
        end
        
        function [dir_index, win] = directionality_index(M, win)
            % dir_index = directionality_index(M, win=[])
            % Wrapper to get directionality index from F.smoothed_direction
            if nargin < 2, win = []; end
            [~, dir_index, win] = M.smoothed_direction(win);
        end
        
        function [dir_sm, dir_index, win, thresh] = smoothed_direction(M, win, thresh)
            % Returns the directions smoothed over <win> second windows
            % where directionality is at least <thresh> (i.e. length of
            % summed vectors is at least <thresh>)
            if nargin < 2 || isempty(win), win = 1; end
            if nargin < 3 || isempty(thresh), thresh = 0.7; end
            
            if median(diff(M.time) < .05)
                dir = exp(1j * M.discharge_directions);
            else
                dir = exp(1j*M.Direction);
            end
            dir_sm = movmean(dir, win, 'omitnan', 'SamplePoints', M.time);
            dir_index = abs(dir_sm);
            N = movsum(isfinite(dir), win, 'samplepoints', M.time);
            dir_index(N == 1) = nan;  % enforce at least 2 discharges in a window to compute the index
            dir_sm(isnan(dir) | dir_index < thresh) = nan;
            dir_sm = angle(dir_sm);
        end
        
        function out = dtheta_dt(~, ~, ~) %#ok<STOUT>
            error('Fix caller: call ''dtheta'' instead of ''dtheta_dt''.');
        end
        function diffs = dtheta(M, win, thresh)
            % Smoothes directions over win and then computes dtheta between
            % consecutive samples 
            if nargin < 2, win = []; end
            if nargin < 3, thresh = []; end
            
            dir = M.smoothed_direction(win, thresh);
            diffs = angle(exp(1j*diff(fillmissing(dir, 'previous'))));
            diffs(isnan(dir(2:end))) = nan;
        end
        
		%% Plotting
        
        function sc = direction_scatter(fit, ax, varargin)
            if nargin < 2; ax = gca; end
            if ~isa(ax, 'matlab.graphics.axis.Axes'), varargin = [ax varargin]; ax = gca; end
            
            [iw_info, mw] = BVNY.get_iw_info(fit.Name);
            if isempty(iw_info), t_offset = 0; else, t_offset = iw_info{mw}.center; end
            
            tt = fit.time - t_offset; tt = [tt; tt; tt];
            yy = fit.Direction; yy = yy + [-1 0 1]*2*pi;
            sc = scatter(ax, tt, yy(:), [], 'filled', ...
                'markerfacealpha', .5, 'tag', 'directions', varargin{:});
            
            % Show the IW
            if ~isempty(iw_info)
                iw_info = iw_info{mw}; 

                sts = ax.NextPlot;
                ax.NextPlot = 'add';
                ylims = 3*[-pi pi];

%                 ln = xline(ax, iw_info.center);  % IW center
                ln = xline(ax, 0);  % IW center

                face_col = .5 * [1 1 1];  % IW bounds
                [v1, v2] = ndgrid(iw_info.range - t_offset, ylims);
                pp = patch(ax, 'vertices', [v1(:) v2(:)], 'faces', [1 2 4 3], ...
                    'facecolor', face_col, ...
                    'facealpha', .3, ...
                    'linestyle', 'none');

                sz = get(sc, 'sizedata');
                dir = scatter(ax, 0, iw_info.phi - unique(fit.RotateBy), ... % IW direction
                    [], 'filled', 'sizedata', max(15, 8*sz(1)), ...
                    'tag', 'iw_dot', 'markerfacealpha', .8);

                ax.NextPlot = sts;
            end
            
            % prettify 
            title(ax, strrep(fit.Name, '_', ' '));
            set(ax, ... 
                'ticklength', [0 0], ... 
                'xgrid', 'on', 'linewidth', 1); 
            set(ax, 'ytick', (-3*pi:pi/2:3*pi));
            lbl = compose('%d\\pi', round(ax.YTick/pi));
            lbl(2:2:end) = {''};
            lbl = strrep(strrep(lbl, '1', ''), '0\pi', '0'); 
            set(ax, 'yticklabel', lbl, ...
                'ygrid', 'on', 'ylim', 1.1*[-pi pi]);
            
        end
        
        function [ax, out] = distance_scatter(M, iw_tpl)
            % iw_tpl = mea.get_IW_templates
            [rho, rho_p, rho_time] = M.correlation(iw_tpl.template);
            mask_ = isfinite(rho_p(:, 1)); 
            scatter(rho_time(mask_), ...
                rho(mask_), ...
                abs(log(rho_p(mask_, 1)+eps)), ...
                rho_p(mask_, 2), 'filled')
            colormap(gca, hsv(max(rho_p(mask_, 2))));
            
            ylim([-1 1])
            xlim(quantile(M.time, [0 1]))
            for tt = iw_tpl.time', xline(tt); end   
            colorbar
            ax = gca;
            ax.CLim = quantile(rho_p(mask_, 2), [0 1]) + [-.5 .5];
            
            out = struct('dist_rho', rho, 'dist_p', rho_p, 'time', rho_time);
        end
        
        function [ax, cmap, bin] = direction_raster(M, ax, dir, cmap)
            % .direction_raster(ax=gca, dir=.Direction, cmap=hsv)
            % All inputs are optional, but ax must come first if given
            
            ax_default = gca;
            dir_default = M.Direction;
            cmap_default = hsv;
            
            switch nargin
                case 1  % F.direction_raster;
                    dir = dir_default;
                    ax = ax_default;
                    cmap = cmap_default;
                case 2
                    if ischar(ax), ax = M.(ax); end
                    if isnumeric(ax)  
                        if size(ax, 2) == 1, dir = ax; cmap = cmap_default; % M.direction_raster(dir)
                        else, cmap = ax; dir = dir_default; % M.direction_raster(cmap)
                        end
                        ax = ax_default;
                    else  % F.direction_raster(ax)
                        dir = dir_default;
                        cmap = cmap_default;
                    end
                case 3
                    if ~isa(ax, 'matlab.graphics.axis.Axes')  % M.direction_raster(dir, cmap)
                        dir = ax; cmap = dir; 
                        ax = ax_default;
                    else
                        if isvector(dir)  % M.direction_raster(ax, dir)
                            cmap = cmap_default;
                        else  % M.direction_raster(ax, cmap)
                            cmap = dir;
                            dir = dir_default;
                        end
                    end
                    if ischar(dir), dir = M.(dir); end
                case 4  % M.direction_raster(ax, dir, cmap)
                    if ischar(dir), dir = M.(dir); end
                    assert( ...
                        isa(ax, 'matlab.graphics.axis.Axes') ...
                        && isvector(dir) ...
                        && size(cmap, 2) == 3, ...
                        'Unexpected inputs. Expected <M.direction_raster(ax, dir, cmap)>.');
                otherwise
                    error('Unexpected inputs. Expected <M.direction_raster(ax, dir, cmap)>.')
            end
            
            
            mask_ = isfinite(dir);
            bin = discretize(dir(mask_), linspace(-pi, pi, size(cmap, 1) + 1));
            yy = [-180; 180] .* ones(2, numel(bin));
            tt = M.time(mask_);
            plot(ax, [tt(:) tt(:)]', yy);
            ax.ColorOrder = cmap(bin, :);
            
            hold(ax, 'on')
            scatter(ax, M.time, dir/pi * 180, 6, [0 0 0], 'filled');
%             plot(ax, M.time, dir/pi*180, 'k.'); 
            hold(ax, 'off');
            
            yticks(-180:90:180);
            title(ax, strrep(M.Name, '_', ' '));
            ax.Tag = 'direction_raster';
            
        end
        
        function h = summary_plot(M, varargin)
            h = figure('name', 'summary_plot');
            T = tiledlayout(h, 5, 1);
            T.TileSpacing = 'compact';
            
            k = 1;  % axis counter
            
            
            ax = nexttile(T, k); 
            k = k+1;
            if median(diff(M.time) < .05)  % if samples are less than 50 ms apart, use discharge times
                dir = M.discharge_directions;
            else
                dir = M.Direction;
            end
            [~, cmap] = M.direction_raster(ax, dir);
            
            title('Direction (\theta)')
            xticklabels(ax, []);
            
            
            % add colorwheel
            pos = T.Children(end).Position;
            ll = pos(1) + pos(3) - pos(4)/3;
            bb = pos(2) + pos(4) - pos(4)/6;
            ww = pos(4)/5;
            hh = ww;
            pax = polaraxes(h, 'units', ax.Units, 'position', [ll bb ww hh]);            
            M.colorwheel(cmap, pax);
            
            
            ax = nexttile(T, k); k = k+1;
            M.direction_raster(ax, 'smoothed_direction');
            title('Smoothed \theta')
            xticklabels(ax, []);
            
            
            win = 1;
            thresh = .7;
            [~, dir_index, ~, thresh] = M.smoothed_direction(win, thresh);
            
            
            ax = nexttile(T, k); k = k+1;  % directionality index
            plot(ax, M.time, dir_index, 'k')
            ylim(ax, [0 1]);
            yline(ax, thresh);
            title('Directionality index')
            xticklabels(ax, []);
            
            
            ax = nexttile(T, k); k = k+1;  % dtheta/dt
            diffs = M.dtheta(win, thresh);
            stem(ax, M.time(2:end), rad2deg(diffs), 'k', 'marker', '.');
            mask_ = abs(diffs) > pi/2;
            tt = M.time(2:end);
            hold on; stem(ax, tt(mask_), rad2deg(diffs(mask_)), 'r', 'marker', '.'); hold off;
            yticks(ax, -180:90:180); grid(ax, 'on')
            ylim([-180 180]);
            title('d(\theta)/dt')
            xlabel('Time [s]')
            ax.Tag = 'dtheta/dt';            
            
            
            ax = nexttile(T, k); k = k+1;  % discharge rate
            dir = M.discharge_directions;
            rate = movsum(isfinite(dir), 1, 'omitnan', 'samplepoints', M.time);
            plot(ax, M.time, rate);
            title('Discharge rate');
            ylabel('Hz');
            
            
            
            linkaxes(T.Children, 'x');
            xlim(ax, [M.time(1), M.time(end)])
            set(findobj(h, 'type', 'axes'), ...
                'xgrid', 'on', 'ygrid', 'on', 'box', 'off');
            title(T, [strrep(M.Name, '_', ' ') ' ' class(M)]);
        end
        
        function [ax, mn, ll, ul] = plot_dirs_with_ci(M, ax, win, units, varargin)
            % Plots directions with 95%CI using a <win> second sliding
            % window.
            % [ax, mn, ll, ul] = plot_dirs_with_ci(M, ax=gca, win=10, units="deg" ::plot args::)
            
            if nargin < 2 || isempty(ax), ax = gca; end
            if nargin < 3 || isempty(win), win = 10; end
            if nargin < 4 || isempty(units), units = "deg"; end
            
            if ~isa(ax, 'matlab.graphics.axis.Axes')  % bump everything forward
                varargin = [units varargin];
                units = win;
                win = ax;
                ax = gca;
            end
            
            if ~isnumeric(win)  % if win isn't numeric, assume it's omitted and bump again
                varargin = [units varargin];
                units = win;
                win = 10;
            end
            
            try % see if units actually units or part of varargin
                units = validatestring(units, ["rad" "deg"]);
            catch ME
                if isempty(varargin), rethrow(ME);
                else
                    varargin = [units varargin];
                    units = "deg";
                end
            end
            
            sts = ax.NextPlot;
            if strcmpi(sts, 'replace'), cla(ax); end  % the patches won't clear the axis first
            
            
            if units == "rad", convert_units = @(x) x; cycle = 2*pi; 
            else, convert_units = @rad2deg; cycle = 360;
            end
            
            switch class(M)
                case 'Delays'
                    mn = M.Direction;
                    ci = M.direction_ci;
                    tt = M.time;
                    
                    % Unwrap so the plot looks nicer and get the right
                    % units
                    ci = convert_units(angle(exp(1j*(ci - mn))));
                    mn = convert_units(unwrap(mn));
                    ci = ci + mn;
                    ul = ci(:, 1); ll = ci(:, 2); 
                    
                otherwise
                    tt = ( M.time(1):.1:M.time(end) )';  % estimate mean and ci at 100 ms intervals
                    [mn, ci] = M.mov_mean_dir(win, tt);
                    mn(isnan(ci) & mn == 0) = nan;
                    
                    % Fill missing ci values, unwrap the directions, and convert to
                    % the requested units.
                    ci = fillmissing(ci, 'constant', pi, 'endvalues', 'none');  % we have no info where ci are nan, so the CI is pi
                    ci = convert_units(ci);
                    mn = convert_units(unwrap(mn));
                    ll = mn - ci;
                    ul = mn + ci;
            end
            
            
                        
            % Plot the mean and 95%CI
            patch_args = varargin;
            patch_args(1:2:end) = strrep(lower(patch_args(1:2:end)), 'color', 'facecolor');
            
            for offset = [-cycle 0 cycle]
                mask_ = isfinite(ll + ul);
                yy = [ll(mask_); flipud(ul(mask_))] + offset;
                xx = [tt(mask_); flipud(tt(mask_))];
                pp = patch(ax, xx, yy, 1, 'displayname', '95%CI', ...
                    'facealpha', .5, 'linestyle', 'none', 'facecolor', [0 0 0], ...
                    patch_args{:}); %#ok<NASGU>
            end
            
            ax.NextPlot = 'add';
            plot(ax, tt, mn + [-cycle 0 cycle], 'k', ...
                'linewidth', 2, 'displayname', 'Mean', varargin{:});
            
            
            % Prettify
            if units == "deg"
                ylim(ax, [-180 180]);
                lbl = cellfun(@(ll) sprintf('%s°', ll), ax.YTickLabel, 'uni', 0);
                ax.YTickLabel = strrep(lbl, '°°', '°');
            else
                ylim(ax, [-pi pi]);
                intv = pi/3;
                tks = -3*pi:intv:3*pi;
                yticks(tks);
                polar_label(ax, 'y');
            end
            title(strrep(M.Name, '_', ' '))
            xlabel('Time [s]')
            ylabel('Direction')
            ax.NextPlot = sts;
            
            
        end

        
        function [ax, mn, ll, ul] = ZZplot_dirs_with_ci(M, ax, win, varargin) % ZZ'd 5/29/21
            % Plots directions with 95%CI using a <win> second sliding
            % window.
            % [ax, mn, ll, ul] = plot_dirs_with_ci(M, ax=gca, win=10, 'units', "deg", 'data', M.Direction, ::plot args::)
            
            default_units = "deg";
            default_ax = gca;
            default_win = 10;
            
            % Parse inputs
            switch nargin
                case 1
                    ax = default_ax; win = default_win;
                case 2
                    if isa(ax, 'matlab.graphics.axis.Axes'); win = default_win;
                    else, win = ax; ax = default_ax;
                    end
                otherwise
                    if isa(ax, 'matlab.graphics.axis.Axes') 
                        if ischar(win), varargin = [win varargin]; win = default_win; end
                    elseif isnumeric(ax), varargin = [win varargin]; win = ax; ax = default_ax;
                    else, varargin = [ax win varargin]; ax = default_ax; win = default_win;
                    end
            end
            ax_status = ax.NextPlot;
            
            % See if units are specified
            units = false(size(varargin));
            units(1:2:end) = contains(lower(varargin(1:2:end)), 'unit');
            if any(units)
                ind = find(units);
                uu = varargin{ind + 1};
                units = validatestring(uu, ["deg" "rad"]);
                varargin(ind + [0 1]) = [];
            else
                units = default_units;
            end
            convert_units = @rad2deg; cycle = 360;
            if units == "rad", convert_units = @(x) x; cycle = 2*pi; end
            
            
            
            tt = ( M.time(1):.1:M.time(end) )';
            [mn, ci] = M.mov_mean_dir(win, tt);  % estimate mean and ci at 100 ms intervals
            nan_ci = isnan(ci);
            ci = fillmissing(ci, 'constant', pi, 'endvalues', 'none');  % we have no info where ci are nan, so the CI is pi
            mn(nan_ci & mn == 0) = nan;
            
            mn = convert_units(unwrap(mn));
            ci = convert_units(ci);
            ll = mn - ci;
            ul = mn + ci;
            
                        
            % Plot the mean and 95%CI
            patch_args = varargin;
            patch_args(1:2:end) = strrep(lower(patch_args(1:2:end)), 'color', 'facecolor');
            
            for offset = [-cycle 0 cycle]
                mask_ = isfinite(ll + ul);
                yy = [ll(mask_); flipud(ul(mask_))] + offset;
                xx = [tt(mask_); flipud(tt(mask_))];
                pp = patch(ax, xx, yy, 1, 'displayname', '95%CI', ...
                    'facealpha', .5, 'linestyle', 'none', 'facecolor', [0 0 0], ...
                    patch_args{:}); %#ok<NASGU>
            end
            
            ax.NextPlot = 'add';
            plot(ax, tt, mn + [-cycle 0 cycle], 'k', ...
                'linewidth', 2, 'displayname', 'Mean', varargin{:});
            
            
            % Prettify
            if units == "deg"
                ylim(ax, [-180 180]);
                lbl = cellfun(@(ll) sprintf('%s°', ll), ax.YTickLabel, 'uni', 0);
                ax.YTickLabel = strrep(lbl, '°°', '°');
            else
                ylim(ax, [-pi pi]);
                intv = pi/3;
                tks = -3*pi:intv:3*pi;
                yticks(tks);
                polar_label(ax, 'y');
            end
            title(strrep(M.Name, '_', ' '))
            xlabel('Time [s]')
            ylabel('Direction')
            ax.NextPlot = ax_status;
            
            
        end
        
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
        
        function dot_size = get_dot_size(M, ax)
            % Get a good dot size to represent an electrode for a scatter
            % plot.
            uu = ax.Units;
            ax.Units = 'points';
            pos = ax.Position;
            ax.Units = uu;
            dot_size = min(pos(3:4))/max(M.GridSize);
            dot_size = dot_size^2;
        end
        
        function ax = plot_residuals(M, varargin)
            [ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
            set(ax.Parent, 'name', 'residuals');
            
            ind = M.time2inds(tt);
            data = squeeze(M.Data(ind, :, :));
            
            % fit a plane
            xx = M.Position(:, 1); yy = M.Position(:, 2);
            zfit = M.Beta(ind, 1) + M.Beta(ind, 2) * xx + M.Beta(ind, 3) * yy;
            
            % Draw a scatter of the data
            dot_size = M.get_dot_size(ax);
			scatter(ax, xx, yy, dot_size, data(:) - zfit, 'filled', 's'); 
            ax.CLim = [-1 1] * quantile(abs(data(:) - zfit), .95);
            colormap(ax, make_diverging_colormap);
            
			title(ax, ["Residuals" ...
                sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', ...
                tt, M.p(ind), M.Magnitude(ind), 'mm/s')]); 
			hold(ax, 'off');
            cb = colorbar;
            title(cb, 'Resid');
            
            xlim(ax, [0 M.GridSize(1) + 1]);
            ylim(ax, [0 M.GridSize(2) + 1]);
			ax.Tag = checkname(['figs' filesep M.Name '_' class(M) '_' ...
                num2str(tt) '_plot_residuals']);
			
        end
        
		function ax = plot3D(M, varargin)
			% ax = plot(obj, t0, ax); % if object has only one time point, t0=obj.time
			[ax, tt] = WaveProp.parse_plot_inputs(varargin{:});
            M_temp = M.sub(M.time2inds(tt));
            
            % Get the data and x-/y-axis positions
            data = squeeze(M_temp.Data);
			N = size(data);
			[p1, p2] = ind2sub(N, find(isfinite(data)));
            
            % Draw a scatter of the data
			scatter3(ax, p1, p2, M_temp.Data(isfinite(data)), [], M_temp.Data(isfinite(data)), 'filled'); 
            
            % Show the fitted plane
			hold(ax, 'on');
			[xx, yy] = ndgrid(1:N(1), 1:N(2));
			zfit = M_temp.Beta(1) + M_temp.Beta(2) * xx + M_temp.Beta(3) * yy;
			im = surf(ax, xx, yy, min(M_temp.Data(:)) * ones(N), data, 'linestyle', 'none'); 
			im.Tag = 'Data';
			surfax = surf(ax, xx, yy, zfit, 'linestyle', 'none', 'facealpha', .5);
			surfax.Tag = 'Fit';
			title(ax, sprintf('T=%.2f\np=%.4g\nspeed=+%.2f%s', M_temp.time, M_temp.p, M_temp.Magnitude, 'mm/s')); 
			hold(ax, 'off');
			ax.Tag = checkname(['figs' filesep M_temp.Name '_' class(M_temp) '_' num2str(tt) '_plot3D']);
			
        end
        
		function [ax, tt] = plot2D(s, varargin)
			% ax = plot2D(obj, tt, ax)
			[ax, tt, directive] = WaveProp.parse_plot_inputs(varargin{:});
            if ndims(s.Data) == 3
				[~, which_t] = min(abs(tt - s.time));
                tt = s.time(which_t);
				data = squeeze(s.Data(which_t, :, :));
				vx = s.Vx(which_t);
				vy = s.Vy(which_t);
                
                pval = s.p(which_t);
                
            else
				data = s.Data;
				vx = s.Vx;
				vy = s.Vy;
                
                pval = s.p;
            end
            phi = atan2(vy, vx);
            speed = norm([vx vy]) * s.MMpElectrode;
            
            % Create an alternate direction using only the first 20% of
            % observed times (gives an idea of how linear the wave
            % trajectory is; if the arrows line up, the trajectory is
            % linear; otherwise, this is a good indication that the waves
            % are not well described by a plane
%             d2 = data;
%             d2(data > quantile(data, .2, 'all')) = nan;
%             [fy, fx] = gradient(d2);
%             Z_ = complex(fx, fy);
%             Z_ = Z_./nanmean(abs(Z_), 'all');
%             v_alt = nanmean(Z_, 'all');
            
            axis(ax, 'square');
            % Get dot_size
            dot_size = s.get_dot_size(ax);
            
            [xx, yy] = find(isfinite(data));
            val = data(isfinite(data)) * 1e3;  % convert from s to ms
            units = 'ms';
            
            % Show the time data
			sc = scatter(ax, xx, yy, dot_size, val, ...
                'filled', 's', directive{:});
            set(ax, 'clim', quantile(sc.CData, [.05, .95]));
			set(ax, 'xtick', [], 'ytick', []);
            
            % Add the arrow
            hold(ax, 'on')
            grid_size = size(data);
            sp = sqrt(vx^2 + vy^2) / (.8*max(grid_size)); % Make the arrow 8 units long
            vx = vx/sp; vy = vy/sp;
            grid_center = grid_size/2 + 1;
            X = grid_center - [vx, vy]/2;  % center the arrow on the array
            qq = quiver(ax, X(1), X(2), vx, vy, s.scale_quiver, ...
                'LineWidth', 2, ...
                'MaxHeadSize', .5, 'color', [0 0 0]);
            
            % Add the alternate arrow (temporary)
%             v_alt = v_alt*.8*max(grid_size);
%             vx = real(v_alt); vy = imag(v_alt);
%             X = grid_center - [vx, vy]/2;  % center the arrow 
%             qq2 = quiver(ax, X(1), X(2), vx, vy, s.scale_quiver, ...
%                 'LineWidth', 2, ...
%                 'MaxHeadSize', .5, 'color', [0 1 0]);
            
            
            hold(ax, 'off')
            ax.Tag = checkname(sprintf("figs/%s_%s_%0.3f_plot2D", s.Name, class(s), tt));
% 			ax.Tag = checkname(['figs' filesep s.Name '_' class(s) '_' num2str(tt) '_plot2D']);
            xlim([0 grid_size(1) + 1]);
            ylim([0 grid_size(2) + 1]);
            axis(ax, 'equal');
            cb = colorbar;
            title(cb, sprintf('TOA [%s]', units))
            title(ax, sprintf('t = %0.3f s', tt));
            lgd_string = ...
                sprintf('\\phi=%0.0f\\circ\n sp=%0.1f mm/s\n p=%0.4g', rad2deg(phi), speed, pval);
            legend(qq, lgd_string)
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
			temp = nan(10);  %max(obj.Position)
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
                ff = validatestring(varargin{ii}, properties(obj));
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
        
        function toa = toa_from_first_N_electrodes(M, N)
            if nargin < 2 || isempty(N), N = 20; end
            
            toa = M.TOA(:, :);
            [~, so] = sort(fillmissing(toa, 'constant', -99), 2, 'descend');  % <fillmissing> because NaNs come first in descent sort
            
            for ii = 1:size(toa, 1)
                toa(ii, so(ii, N+1:end)) = nan;
            end
            
        end
        
		function M = refit_data(M)
			M.Inds = [];  % Remove time resampling
            M.ResultOfFit = [];
            M.ResultOfFit;
% 			data = M.Data;
%             [px, py] = find(squeeze(any(isfinite(data))));
%             M.Position = [px, py];
%             
%             N = size(data, 1);
%             locs = sub2ind(size(data, [2 3]), M.Position(:, 1), M.Position(:, 2));
%             
%             for ii = 1:N
%                 dat = data(ii, locs);
%                 [V, M.p(ii), beta] = M.fit_data(dat, M.Position, M.Quadratic);
%                 M.Vx(ii) = V(1);
%                 M.Vy(ii) = V(2);
%                 M.Beta(ii, :) = beta;		
%                 M.Curvature(ii) = WaveProp.curvature_from_beta(beta);
%                 ZZ = complex(V(1), V(2));
%                 ZZ(abs(ZZ) > 1e6) = nan;
%                 M.Magnitude(ii) = abs(ZZ);
%                 M.Direction(ii) = angle(ZZ);
%                 
%             end
%             warning(W);
            
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
            F = copy(obj);
            F.Inds = ind;
        end
        
		function F = split(obj)
			N = min(numel(obj.time), numel(obj.Vx));
			F(N) = obj.sub(N);
			for ii = 1:N
				F(ii) = obj.sub(ii);
			end
        end
        
		function d = diff(s, other)
            i0 = other.Inds;
			other = other.resample_t0(s.time);
			d = angle(exp(1j*(s.Direction - other.Direction)));
            other.Inds = i0;  % Return inds to what they were (these are handles now. Without this line, this operation changes <other>)
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
    
    methods  % Reload
        function obj = reload(obj, S)
            % I renamed some fields on 5/10/21; this is just to make sure
            % everything gets reloaded correctly. Won't be necessary after
            % the next re-run
            if isstruct(S) && isfield(S, 'Data')
                obj.GridSize = size(S.Data, [2 3]);
                [px, py] = ndgrid(1:obj.GridSize(1), 1:obj.GridSize(2));
                obj.Position = [px(:), py(:)];
                obj.TOA = S.Data(:, :);
            else
                obj = S;
            end
        end
    end
		
	methods (Static)
        function obj = loadobj(S)
            
            if isstruct(S)
                obj = WaveProp;
                obj = reload(obj, S);
            end
            
        end
        
		function [ax, tt, directive] = parse_plot_inputs(ax, tt, varargin)
            directive = {}; 
            deft_tt = 0;
            deft_ax = @() gca;  % axes(figure);  % just made 70 figures because I forgot to give an axis... probably going to crash matlab
            switch nargin
                case 0
                    ax = deft_ax(); tt = deft_tt; % default to new figure
                case 1
                    if isa(ax, 'matlab.graphics.axis.Axes'), tt = deft_tt; 
                    else, tt = ax; ax = deft_ax();
                    end
                otherwise
                    if isa(ax, 'matlab.graphics.axis.Axes') 
                        if isnumeric(tt), directive = varargin;
                        else, directive = [tt, varargin]; tt = deft_tt; 
                        end
                    elseif isnumeric(ax)
                        if isa(tt, 'matlab.graphics.axis.Axes')
                            directive = varargin;
                            tt_ = tt; tt = ax; ax = tt_;
                        else
                            directive = [tt varargin];
                            tt = deft_tt;
                        end
                    else
                        ax = deft_ax(); tt = deft_tt; directive = varargin;
                    end
            end
			
        end
		
        function pax = colorwheel(cmap, pax)
            if nargin < 2, pax = []; end
            if nargin < 1 || isempty(cmap), cmap = hsv; end
            
            if isempty(pax)
                if isa(gca, 'matlab.graphics.axis.PolarAxes'), 
                    pax = gca; 
                else
                    h = figure('name', 'cwheel', 'units', 'inches', ...
                        'position', [0 0 1 1], 'resize', 'off');
                    pax = polaraxes(h);
                end
            end
            
            theta_lims = linspace(-pi, pi, size(cmap, 1));
            ph =@(ii) polarhistogram(pax, theta_lims(ii), theta_lims(ii:ii+1), ...
                'linestyle', 'none', 'facecolor', cmap(ii, :), 'facealpha', 1);
            pax.NextPlot = 'add';
            arrayfun(@(ii) ph(ii), 1:size(cmap, 1) - 1);
            pax.NextPlot = 'replace';
            set(pax, 'rtick', [], 'rlim', [0 2], 'nextplot', 'replace', ...
                'thetatick', -pi/2:pi/2:pi, 'thetalim', [-pi pi], ...
                'thetaaxisunits', 'degrees', 'color', 'none', ...
                'rdir', 'reverse');
            pax.ThetaTickLabel = cellfun(@(s) [s '\circ'], pax.ThetaTickLabel, 'uni', 0);
            
            pax.Tag = 'colorwheel';

            
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
        
        function [V, p, beta] = fit_data(data, position, quadratic)
            % At one point I was messing with fitting quadratics. This will
            % probably never happen again, but who knows.
            if nargin < 3, quadratic = false; end
            if quadratic, [V, p, beta] = WaveProp.fit_quadratic(data, position);
            else, [V, p, beta] = WaveProp.fit_plane(data, position);
            end
        end
        
        function [V, p, beta] = fit_quadratic(data, position)
            % Fit a quadratic plane to the data. Return nan if not enough
            % data or if beta are zero. 
            
			V = nan(1, 2); p = nan; beta = nan(1, 6);
            if numfinite(data(:)) < 3, return; end
            
            X = [position prod(position, 2) position.^2];
            
            % fit the delay vs two-dimensional positions
            [beta,stats] = robustfit(X, data, 'fair'); 
            
            % perform F test that last two coefs are both 0.
            H = [0 1 0; 0 0 1];  % These are defaults, but here for clarity
            c = [0 ; 0];
            p = linhyptest(beta, stats.covb, c, H, stats.dfe);  
			
			% beta is invalid if nan or if the slope is 0 in both directions
			invalid = all(beta(2:3).^2 < eps);
			if invalid, return; end		
            
            V = pinv(beta(2:3));
        end
        
		function [V, p, beta, stats] = fit_plane(data, position)
			% Fit a 2D plane to the data. Return nan if not enough data or
			% beta are zero.
			V = nan(1, 2); p = nan; beta = nan(1, 3); stats = nan;
            if numfinite(data(:)) < 4, return; end
            
            % fit the delay vs two-dimensional positions
            W = warning;
            warning('off', 'stats:statrobustfit:IterationLimit');
            warning('off', 'stats:robustfit:RankDeficient');
            warning('off', 'MATLAB:singularMatrix');
            [beta, stats] = robustfit(position, data, 'fair'); 
            
            % perform F test that last two coefs are both 0.
            H = [0 1 0; 0 0 1];  % These are defaults, but here for clarity
            c = [0 ; 0];
            p = linhyptest(beta, stats.covb, c, H, stats.dfe);  
            
			warning(W);
            
			% beta is invalid if nan or if the slope is 0 in both directions
			invalid = all(beta(2:3).^2 < eps);
			if invalid, return; end			

			V = pinv(beta(2:3));
        end
        
		function [reduced_data, G] = use_largest_cluster(data, position)
			% Use the largest cluster of data points for lfp methods
            % Assumes first dimension is time; if 2D, <position> must be
            % given.
            
            if nargin < 2
                assert(ndims(data) == 3, 'Data must be 3d or electrode position must be given.')
                [px, py] = find(ones(size(data, [2 3])));
                position = [px, py];
            end
            
            reduced_data = nan(size(data));
            G = nan(size(data));
			warning('off', 'stats:clusterdata:MissingDataRemoved');
            
            for ii = 1:size(data, 1)
                temp_dat = data(ii, :);
                X = [position, normalize(temp_dat(:))];
                if numfinite(X(:, end)) < 3, continue; end

                cc = clusterdata(X, ...
                    'distance', 'euclidean', ... % euclidean distance
                    'criterion', 'distance', ...
                    'linkage', 'single', ...   % Nearest neighbor
                    'cutoff', sqrt(2 + .5));  % ... must be connected and values differ by less than .5 sd
                [~, largest] = max(histcounts(cc, max(cc)));
                mask = cc == largest;
%                 reduced_data = data;
                temp_dat(~mask) = nan;
                reduced_data(ii, :) = temp_dat;
                G(ii, :) = cc;
            end
			
        end
        
		function S = resize_obj(s, to_struct)
			if nargin < 2, to_struct = false;  end
			
            if numel(s) == 1; S = s; return; end
             % If all positions match, continue (if not, do not combine
             % these objects)
            temp_position = cat(3, s.Position);
            assert(all(temp_position(:, :, 1) == temp_position(:, :, 1:end), 'all'), ...
                'Field <Position> is not the same for all objects');
            
			S = copy(s(1));
            
            if to_struct
				S = struct;
				for p = string(properties(s(1))')
					S.(p) = []; 
				end
            end
            
            toa_new = cat(1, s.TOA);
            S.TOA = toa_new(:, :);
            S.t0 = cat(1, s.t0);            
            S.Inds = [];
            
            mco = metaclass(s).PropertyList;
            FF = {mco.Name};
            FF(cat(1, mco.Dependent) | cat(1, mco.Hidden) | cat(1, mco.Transient)) = [];
%             out = struct;
            for f = string(FF)
                if ismember(lower(f), ["toa", "t0", "position"]), continue; end
                new_val = cat(1, s.(f));
                siz = size(new_val);
                % Get unique values
                if isvector(new_val), uu = unique(new_val);
                else, uu = unique(new_val(:, :), 'rows');
                end
                
                if size(uu, 1) == 1
                    S.(f) = reshape(uu, [1, siz(2:end)]);
                else
                    S.(f) = new_val;
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
            %          'MG49_Seizure43'), (3) a cell list of files from the
            %          WaveFits folder or (4) it can be empty in which case
            %          fits from all seizures in SeizureInfo will load
            
            if nargin < 2, metrics = []; end
            if nargin < 1 || isempty(mea)
                sz = SeizureInfo;
                if isempty(metrics), metrics = findsharedmetrics_(sz); end
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
            elseif iscell(mea)
                % Load a list of files in WaveFits folder
                f0 = cellfun(@(s) WaveProp.load(s, metrics), mea);
                for ff = fieldnames(f0)'
                    fits.(ff{:}) = cat(1, f0.(ff{:}));
                end
            elseif isa(mea, 'MEA')
                fits = foldercontents2struct_(['WaveFits/' mea.Name], metrics);
                fits.Name = {mea.Name};  % convert to cell to match previous version of load
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