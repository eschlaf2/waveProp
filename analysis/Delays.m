classdef Delays < WaveProp
    
	properties
        CI_delays  % confidence intervals surrounding the delays on each electrode
        CI_dir  % generate conficence intervals for the direction estimates
        
		HalfWin = 5
		FBand = [1 13]
		MinFreq = 3
		Reference string = 'center'
		MaskBy string = 'longest'
		Freq = [nan; nan]
		Type string = 'group'
% 		SamplingRate
		
		BW = 2 % Bandwidth
		NTapers  % # of tapers
		FFTPad = 4
		CohConf = 0.05
		T
        NCh
		
	end

	properties (Transient = true, Access = private)
		window
        
        % Used to get the TOA, but no need to save afterward
        Time  % This will just be indices unless given as input
        SamplingRate = 1  % Assumes 1 sample per second unless otherwise given
        
    end
    
    properties (Dependent = true, Access = public)
		RefTrace
	end
    
    
	methods  % Constructor
		function D = Delays(mea, t0, varargin)	
		% Computes the delays at time t0 for the data in <obj>. 
		% Inputs: 
		%    obj: struct with fields 
		%        Data (required): TxN matrix (T=time, N=channels)
		%        Time (optional): length T vector of time points
		%        SamplingRate (optional): scalar sampling rate of the data
		%            *** one of Time or SamplingRate must be present ***
		%        Position (required): Nx2 matrix with position of channels
		%    t0: time at which to test the delays (centered window)
		%
			
        % Early days messing with classes... Probably better to rewrite
        % from scratch if you get irritated.
        
			if nargin < 1, return, end
			if nargin < 2, t0 = []; end
			
            % If multiple times are given, loop through and then combine
            % objects. This is a really slow way to do this, but things are
            % working and I don't want to reorganize everything now...
            % Tolerate the slow and then optimize later if you keep using
            % this.
            if numel(t0) > 1
                N = numel(t0);
                D0(N) = Delays(mea, t0(N), varargin{:});
                tic
                for ii = 1:N-1
                    if ~mod(ii, 10)
                        fprintf('Computed delay %d/%d (%0.2f seconds)\n', ii-1, N, toc); 
                        tic; 
                    end
                    D0(ii) = Delays(mea, t0(ii), varargin{:});
                end
                
                D = WaveProp.resize_obj(D0);
                return
            end
            
			try signal = mea.Data; catch, signal = mea.data; end
			for ff = string(fieldnames(mea)')
				switch lower(ff)
					case 'time'
						D.Time = mea.(ff);
					case 'name'
						D.Name = mea.(ff);
					case 'position'
						D.Position = mea.(ff);
					case 'samplingrate'
						D.SamplingRate = mea.(ff);
						if isempty(D.Time), D.Time = 1:size(signal, 1)' / mea.(ff); end
                    case 'patient'
                        D.Patient = mea.(ff);
                    case 'seizure'
                        D.Seizure = mea.(ff);
                    case 'gridsize'
                        D.GridSize = mea.(ff);
					otherwise
						continue
				end
			end
			if isempty(t0), t0 = mean(D.Time); end
            
			D.t0 = t0;
            D.NCh = length(D.Position);
			D = D.parse_inputs(varargin{:});			
			
            % Get the data window
            t_inds = abs(D.Time - D.t0) <= D.HalfWin;
			D.T = range(D.Time(t_inds));
			D.window = signal(t_inds, :);
            
            % compute the TOA
            [coh, phi, freq, coh_conf] = ...
				D.compute_coherence(D.window, D.RefTrace, D.params);
			coh(coh < coh_conf) = nan;
			[phi, freq_out] = D.mask_phi(coh, phi, freq);
            [delay, delay_ci] = D.regress_delay(phi, freq, D.Type);
% 			delay = D.phi2delay(phi, freq);

            if isempty(delay)
                delay = nan;
            end
            
            % Put time in the first dimension, electrodes in the second
            % dimension, other stuff in later dimensions
            D.TOA = reshape(-delay, [1, numel(delay)]);
            D.Freq = reshape(freq_out, [1 size(freq_out)]);
            D.CI_delays = reshape(delay_ci, [1 size(delay_ci)]);
	
        end
        

    end
    

	methods  % coherence parameters
		function sr = get.SamplingRate(obj)
			sr = obj.SamplingRate;
			if isempty(sr)
				sr = 1 / diff(obj.Time(1:2));
			end
		end
		
		function mf = get.MinFreq(obj)
			mf = obj.MinFreq;
			if isempty(mf)
				mf = max(obj.DefaultMinFreq, min(3 / obj.HalfWin, range(obj.FBand)));
			end
		end
		
		function ntapers = get.NTapers(obj)
			ntapers = obj.NTapers;
			if isempty(ntapers)
				ntapers = floor(2 * (obj.T .* obj.BW)) - 1;  % one less than the shannon number
			end
            if numel(ntapers) == 1 && ntapers < 5
                ntapers = 5;
                warning('Number of tapers is changed to 5 (from %d; new BW=%0.2f).', ntapers, (ntapers + 1) / 2 / obj.T); 
            end
		end
    end
    
	methods  % getters for dependent properties
		
		function reference = get.RefTrace(D)
			switch D.Reference
				case 'mean'
					reference = Delays.mean_signal(D.window);
				case 'center'
					reference = Delays.central_electrode(D.window, D.Position);
				otherwise
					error('Reference %s not recognized', D.Reference);
			end
        end
        
		function p = params(D)
			% Create the params struct that will be passed to chronux
			% function conherencyc
			p.tapers = [D.T .* D.BW, D.NTapers];  % ... time-bandwidth product and tapers.
			p.Fs = D.SamplingRate; % ... sampling rate
			p.pad = D.FFTPad;                 % ... 2^(ceil(log2(T)) + pad)
			p.fpass = D.FBand;            % ... freq range to pass
			p.err = [1 D.CohConf];          % ... theoretical error bars, p=0.05.
			p.T = D.T;
		end

	end
	
	properties (Hidden = true)
		DefaultMinFreq = 3
		ParamNames = ["Position" "HalfWin" "FBand" "MinFreq" "Reference" ...
			"MaskBy" "BW" "NTapers" "FFTPad" "CohConf" "Type"]
	end
	
	methods
		
        function ci = bootstrap_direction_ci(D, Nboot, correction)
            % Resamples from the CI of the delays on each channel and then
            % recomputes direction based on the resampled delays.
            %
            % This is not a great method, but it's interesting. The CI end
            % up biased, though, because (I assume) we don't account for the
            % covariance structure between in the electrodes. I "fixed it"
            % by subtracting the mean direction from the bootstrap samples
            % and then re-centering around the reported direction. This is
            % probably not really what we want to look at which is why I
            % say it's not a great method. 
            
            if nargin < 2 || isempty(Nboot), Nboot = 1000; end
            if nargin < 3 || isempty(correction), correction = true; end
            
            toa_est = D.TOA;
            delay_ci = D.CI_delays;
            N = size(delay_ci, 1);
            Nch = size(toa_est, 2);
            dir_est = D.Direction;
            
            ci = nan(N, 2);
            for ii = 1:N
                if isnan(dir_est(ii)), continue; end
                
                new_dir = nan(Nboot, 1);
                for jj = 1:Nboot
                    
                    % Get a set of possible TOAs based on CI for each
                    % channel. 
                    toa_sample = ...
                        rand(1, Nch) .* diff(delay_ci(ii, :, :), [], 3) ...
                        - delay_ci(ii, :, 1);
                    [V, p_, beta] = D.fit_plane(toa_sample, D.Position);
                    new_dir(jj) = atan2(V(2), V(1));
                end
                
                % center directions around 0 so quantiles are estimated
                % correctly.
                dir0 = circ_mean(new_dir);
                new_dir = angle(exp(1j* (new_dir - dir0)));  
                
                % Without this, ci are not centered around mean. I think is
                % because bootstrapping using the 95%CI of the delays
                % doesn't account for covariance between the electrodes.
                if correction, dir0 = dir_est(ii); end
                
                % Compute CI and return to proper baseline
                ci(ii, :) = quantile( new_dir, [.025, .975]) + dir0;
                ci(ii, :) = fix_angle(ci(ii, :));  % keep in [-π,π]
            end
            
        end
        
		function [phi, freq_out] = mask_phi(D, coh, phi, freq)
			switch lower(D.MaskBy)
				case 'longest'
					mask = D.reduce_data(isfinite(coh), D.MinFreq / diff(freq(1:2)));					
				case 'highest'
					mask = D.coh_mask(coh, D.MinFreq / diff(freq(1:2)));
				otherwise 
					error('Value of ''MaskBy'' unrecognized.')
			end
			phi(~mask) = nan;
			freq = repmat(freq(:), 1, size(phi, 2));
			freq(~mask) = nan;
			freq_out = [min(freq); max(freq)]';  % return electrodes along the first dimension
        end
        
		function [delay, delay_ci] = phi2delay(D, phi, freq)
			switch lower(D.Type)
				case 'group'
                    [delay, delay_ci] = D.regress_delay(phi, freq);
% 					delay = D.group_delay(phi, freq);  % compute delays on each electrode based on coherence
				case 'phase'
					delay = D.phase_delay(phi, freq);
                case 'regress'
                    delay = D.regress_delay(phi, freq);
				otherwise
					error("'Type' not recognized. Must be 'group' or 'phase'.")
			end
        end
        
        function obj = reload(obj, S)
            % Made some changes to the WaveProp properties. This is to make
            % sure that there are no errors in loading after the changes. 
            obj = reload@WaveProp(obj, S);
        end
	end
	
	
	%% Static methods
	methods (Static)
		
        
        
        function obj = loadobj(S)
            if isstruct(S)
                obj = Delays;
                for ff = string(fieldnames(S)')
                    if ismember(ff, {'ParamNames', 'WPParamNames'}), continue; end
                    obj.(ff) = S.(ff); 
                end
                obj = reload(obj, S);
            else
                obj = S;
            end
        end
        
        function [delay, delay_ci] = regress_delay(phi, freq, delay_type)
            % computes group delay be fitting a line to phi(f) and
            % returning slope

            if nargin < 3, delay_type = 'group'; end
            delay_type = validatestring(delay_type, ["group" "phase"]);
            
            
            [d1, dint, ~, ~, stats] = arrayfun(@(ii) ...
                regress(phi(:, ii), [ones(size(freq)); freq]'), ...
                1:size(phi, 2), 'uni', 0);
            d1 = cat(2, d1{:})';
            dint = permute(cat(3, dint{:}), [3 2 1]);  % put electrode dimension first; then bounds, then variable (intercept and slope). dint(trode, 1, 2) is the lower bound of the slope
            stats = cat(1, stats{:});
            
            if delay_type == "group", D = [d1(:, 2) dint(:, :, 2)];  % use slope to estimate delay
            elseif delay_type == "phase", D = [d1(:, 1) dint(:, :, 1)];  % use intercept to estimate delay
            end
            
            pval = stats(:, 3);
            
            D(pval >= .05, :) = nan;
            delay = -D(:, 1) ./ (2*pi);
            delay_ci = squeeze(-D(:, 2:3) ./ (2*pi));
        end
        
		function delay = group_delay(phi, freq)
             % computes group delay by taking mean dphi/df
            df = diff(freq);
			delay = nanmean(-diff_phase(phi) ./ df(:)) / (2*pi);
        end
        
		function delay = phase_delay(phi, freq)
			delay = nanmean(-phi ./ freq(:)) / (2*pi);
        end
        
		function mask = coh_mask(coh, min_length)
			finite_coh = isfinite(coh);
			N = size(coh, 2);
			groups = cumsum(~finite_coh);
			groups(~finite_coh) = nan;
			
			mask = false(size(coh));
			for ii = 1:N
				G = findgroups(groups(:, ii));
				mn = splitapply(@mean, coh(:, ii), G);
				nml = splitapply(@numel, coh(:, ii), G);
				[~, so] = sort(mn, 'descend');
				too_short = nml(so) < min_length;
				so(too_short) = [];
				if isempty(so), continue, end
				mask(:, ii) = G == so(1);
			end
        end
        
		function mask = reduce_data(valid, min_length)
% 			valid = isfinite(valid);
			N = size(valid, 2);
% 			data_r = nan(size(data));
			mask = false(size(valid));
			for ii = 1:N
				starts = find(diff([0; valid(:, ii)]) > 0);
				ends = find(diff([valid(:, ii); 0]) < 0);
				[length, ind] = max(ends - starts);
				if length+1 < min_length, continue, end
				mask(starts(ind):ends(ind), ii) = true;
% 				data_r(starts(ind):ends(ind), ii) = ...
% 					data(starts(ind):ends(ind), ii);
			end
        end
        
		function [coh, phi, freq, coh_conf] = compute_coherence(window, reference, params)

			window = window - nanmean(window);
			N = size(window, 2);
			reference = reference(:) .* ones(1, N);
			[coh, phi, ~, ~, ~, freq, coh_conf, ~] = ...
				coherencyc(window, reference, params);
			
		end
		
		function reference = mean_signal(window)
			window = window - nanmean(window);
			reference = nanmean(window, 2);
			
        end
        
		function reference = central_electrode(window, position)
			center = mean(position);
			dist2center = sum((position - center).^2, 2);
			[~, centralI] = min(dist2center);
			reference = window(:, centralI);
		end

	end
	
	
	
end
