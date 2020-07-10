classdef Delays < WaveProp
	methods  % Constructor
		function D = Delays(obj, t0, varargin)	
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
			
			if nargin < 1, return, end
			if nargin < 2, t0 = []; end
			
			try data = obj.Data; catch, data = obj.data; end
			for ff = string(fieldnames(obj)')
				switch lower(ff)
					case 'time'
						D.Time = obj.(ff);
					case 'name'
						D.Name = obj.(ff);
					case 'position'
						D.Position = obj.(ff);
					case 'samplingrate'
						D.SamplingRate = obj.(ff);
						if isempty(D.Time)
							D.Time = 1:size(data, 1)' / obj.(ff);
						end
					otherwise
						continue
				end
			end
			if isempty(t0), t0 = mean(D.Time); end
			D.t0 = t0;
			D = D.parse_inputs(varargin{:});
			[D.window, D.T] = D.get_window(data);
			[data, D.Freq] = D.get_data();
			if numel(unique(data)) < 3, return, end
			D.Data = -data;
			D = D.compile_results(); 	
		end

	end
	properties
		HalfWin = .5
		FBand = [1 13]
		MinFreq = []
		Reference = 'center'
		MaskBy = 'longest'
		Freq = [nan; nan]
		Type = 'group'
		SamplingRate
		
		BW  % Bandwidth
		NTapers  % # of tapers
		FFTPad = 4
		CohConf = 0.05
		T
		
	end

	properties (Transient = true, Access = private)
		window
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
		function bw = get.BW(obj)
			bw = obj.BW;
			if isempty(bw)
				switch obj.HalfWin
					case 5
						bw = 2;
					case 0.5
						bw = 3;
					otherwise
						bw = 3 / obj.T;  % This will generate 5 tapers
				end
			end
		end
		
		function ntapers = get.NTapers(obj)
			ntapers = obj.NTapers;
			if isempty(ntapers)
				ntapers = floor(2 * (obj.T .* obj.BW)) - 1;  % one less than the shannon number
			end
			if numel(ntapers) == 1 && ntapers < 5, warning('Number of tapers is less than 5: %d', ntapers); end
		end
	end
	properties (Dependent = true, Access = private)
		RefTrace
		NCh
	end
	methods  % getters for dependent properties
		function N = get.NCh(D)
			N = size(D.Data);
			N = prod(N(2:end));
		end
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
			"MaskBy" "BW" "NTapers" "FFTPad" "CohConf"]
	end
	
	methods
		
        function freq = get.Freq(D)
            freq = D.Freq;
            switch class(freq)
                case 'cell'
                    if numel(freq) == numel(D.time)
                        empties = cellfun(@isempty, freq);
                        fill_with = nan(size(freq{find(~empties, 1)}));
                        freq(empties) = {fill_with};
                        dim = find(size(fill_with) == 1, 1);
                        freq = cat(dim, freq{:});
                        if dim ~= 1, 
                            freq = shiftdim(freq, dim - 1); 
                        end
                    end
            end
        end
		function [data_win, T] = get_window(D, data)
			t_inds = abs(D.Time - D.t0) <= D.HalfWin;
			T = range(D.Time(t_inds));
			data_win = data(t_inds, :);
		end

		function [delay, freq_out] = get_data(D)
			
			[coh, phi, freq, coh_conf] = ...
				D.compute_coherence(D.window, D.RefTrace, D.params);
			coh(coh < coh_conf) = nan;
			[phi, freq_out] = D.mask_phi(coh, phi, freq);
			delay = D.phi2delay(phi, freq);
			
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
			freq = repmat(freq(:), 1, D.NCh);
			freq(~mask) = nan;
			freq_out = [min(freq); max(freq)];
		end
		function delay = phi2delay(D, phi, freq)
			switch lower(D.Type)
				case 'group'
					delay = D.group_delay(phi, freq);      % compute delays on each electrode based on coherence
				case 'phase'
					delay = D.phase_delay(phi, freq);
				otherwise
					error("'Type' not recognized. Must be 'group' or 'phase'.")
			end
		end
	end
	
	
	%% Static methods
	methods (Static)
		
		function delay = group_delay(phi, freq)
			df = diff(freq);
			delay = nanmean(-diff_phase(phi) ./ df(:));
		end
		function delay = phase_delay(phi, freq)
			delay = nanmean(-phi ./ freq(:));
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
% 		function [params, compute_inds] = set_coherence_params(Time, T, band)
% 
% 		% 	band = [1 13];                  % Select a frequency range to analyze
% 			W = 3/T;                          % Bandwidth
% 			NTAPERS = round(2*(T * W) - 1);        % Choose the # of tapers.
% 			OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)
% 
% 			samplingRate = 1 / diff(Time(1:2));
% 
% 			params.tapers = [T * W, NTAPERS];  % ... time-bandwidth product and tapers.
% 			params.Fs = samplingRate; % ... sampling rate
% 			params.pad = 4;                 % ... 2^(ceil(log2(T)) + pad)
% 			params.fpass = band;            % ... freq range to pass
% 			params.err = [1 0.05];          % ... theoretical error bars, p=0.05.
% 			params.T = T;
% 
% 			nsamp = numel(Time);                                % Total number of samples
% 			step = OVERLAP_COMPLEMENT * samplingRate;           % Compute coherence every OVERLAP_COMPLEMENT
% 			lastSamplePoint = nsamp - ceil(T * samplingRate) + 1;     % Leave a long enough window at the end to calculate coherence
% 			compute_inds = round(1 : step : lastSamplePoint);    
% 
% 		end
