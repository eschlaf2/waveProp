classdef Delays < WaveProp
	
	properties
		HalfWin = .5
		FBand = [1 13]
		MinFreq = []
		Reference = 'center'
		MaskBy = 'longest'
		Freq
		Type = 'group'
	end
	
	properties (Hidden = true)
		DefaultMinFreq = 3
		ParamNames = ["Position" "HalfWin" "FBand" "MinFreq" "Reference" "MaskBy"]
	end
	
	methods
		
% 		function obj = Delays(mea, t0, position, T, half_win, fband, min_freq)
		function obj = Delays(mea, t0, varargin)			
			
			if nargin < 1, return, end
			if isa(mea, 'MEA')
				obj.Time = mea.Time;
				obj.Name = mea.Name;
				obj.t0 = t0;
				obj.HalfWin = mea.params.T / 2;
				obj.FBand = mea.params.delay_band;

				obj.Position = mea.Position;
				obj = obj.parse_inputs(varargin{:});
				[window, T] = obj.get_window(mea);
			else
				obj = obj.parse_inputs(varargin{:});
				window = mea;
			end
			if isempty(obj.MinFreq)
				obj.MinFreq = max(obj.DefaultMinFreq, min(3 / obj.HalfWin, range(obj.FBand)));
			end
			
			
			[data, obj.Freq] = obj.get_data(window, T);
			if numel(unique(data)) < 3, return, end
			obj.Data = data;
			obj = obj.compile_results; 	
		end
		
		function [window, T] = get_window(obj, mea)
			t_inds = abs(mea.Time - obj.t0) <= obj.HalfWin;
			T = obj.Time(t_inds);
			data = mea.Data;
			window = data(t_inds, :);
		end
		function [data, freq_out] = get_data(s, window, T)
			
% 			T = s.Time(s.t_inds);
% 			data = mea.Data;
% 			data = data(s.t_inds, :);
% 			[window, T] = s.get_window(mea);
% 			time_series_data = mea.filter( ...
% 				data(s.t_inds, :), ...
% 				mea.SamplingRate, mea.params.delay_band);
% 			D = [];
			[params, ~] = s.set_coherence_params( ...
				T, range(T), s.FBand);
% 			[~, center] = min(sum((s.Position - mean(s.Position)).^2, 2));         % find the most central electrode
			switch s.Reference
				case 'mean'
					reference = Delays.mean_signal(window);
				case 'center'
					reference = Delays.central_electrode(window, s.Position);
				otherwise
					error('Reference %s not recognized', s.Reference);
			end
% 			[coh, phi, freq, coh_conf] = ...
% 				compute_coherence(window, params, 'pairs', center);          % compute the coherence over the selected interval
			[coh, phi, freq, coh_conf] = ...
				s.compute_coherence(window, reference, params);
			coh(coh < coh_conf) = nan;
			switch lower(s.MaskBy)
				case 'longest'
					mask = s.reduce_data(isfinite(coh), s.MinFreq / diff(freq(1:2)));					
				case 'highest'
					mask = s.coh_mask(coh, s.MinFreq / diff(freq(1:2)));
				otherwise 
					error('Value of ''MaskBy'' unrecognized.')
			end
			phi(~mask) = nan;
			switch lower(s.Type)
				case 'group'
					data = s.group_delay(phi, freq);      % compute delays on each electrode based on coherence
				case 'phase'
					data = s.phase_delay(phi, freq);
				otherwise
					error("'Type' not recognized. Must be 'group' or 'phase'.")
			end
			freq = repmat(freq(:), 1, size(window, 2));
			freq(~mask) = nan;
			freq_out = [min(freq); max(freq)];
		end
	end
	
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
		function [params, compute_inds] = set_coherence_params(Time, T, band)

		% 	band = [1 13];                  % Select a frequency range to analyze
			W = 3/T;                          % Bandwidth
			NTAPERS = round(2*(T * W) - 1);        % Choose the # of tapers.
			OVERLAP_COMPLEMENT = 1;         % T - OVERLAP (s)

			samplingRate = 1 / mean(diff(Time));

			params.tapers = [T * W, NTAPERS];  % ... time-bandwidth product and tapers.
			params.Fs = samplingRate; % ... sampling rate
			params.pad = 4;                 % ... 2^(ceil(log2(T)) + pad)
			params.fpass = band;            % ... freq range to pass
			params.err = [1 0.05];          % ... theoretical error bars, p=0.05.
			params.T = T;

			nsamp = numel(Time);                                % Total number of samples
			step = OVERLAP_COMPLEMENT * samplingRate;           % Compute coherence every OVERLAP_COMPLEMENT
			lastSamplePoint = nsamp - ceil(T * samplingRate) + 1;     % Leave a long enough window at the end to calculate coherence
			compute_inds = round(1 : step : lastSamplePoint);    

		end

	end
	
	
	
end