classdef (HandleCompatible) MEA < matlab.mixin.Heterogeneous & handle

	properties 
		BadChannels
		Data
		Duration
		Map
		Name
		Padding
		Path = 'MG49/MG49_Seizure43_Neuroport_10_10.mat'
		Position
		Time
		SamplingRate = 1000
		event_inds
	end
	
	properties (Hidden = true)
		locs
		patient
		seizure
		Raw
		AllTime
		MaxDescentData
		PairwiseCoherence
	end
	
	properties (Access = public, Hidden = true)
		params = init_mea_params();
		SRO
		lfp
		mua
		artefacts
		lfp_lo
		custom
		firing_rate
		mua_events
		wave_fits
		event_times
		metrics
% 		flow
% 		metrics
		
		Spectrum
		BSI
		WaveTimes
		ExcludedChannels
		Fits
		Coherent
		PeakRatio
	end
	
	properties (Access = public, Dependent = true)
		skipfactor
		Center
		Active
	end
	

	methods % init, setters and getters
		
		function mea = MEA(path, SR)
			
			if exist('SR', 'var'), mea.SamplingRate = SR; end
			if exist('path', 'var'), mea.Path = path; end
			mea.load;
		end
		
		function mea = load(mea)
			temp = load(mea.Path);
			mea.SRO = temp.SamplingRate;
			mea.AllTime = temp.Time;
			mea.Raw = temp.Data;

			for f = fieldnames(temp)'
				switch f{:}
					case {'Data', 'SamplingRate', 'Time', 'Path'}
						continue
					case fieldnames(mea)
						mea.(f{:}) = temp.(f{:});
					otherwise
						continue
				end
			end
			
			[mea.patient, mea.seizure] = mea.get_info(mea.Path);
			
		end
		function name = get.Name(mea)
			name = [mea.patient '_Seizure' mea.seizure];
		end
		function active = get.Active(mea)
			active = (mea.Time > 0) & (mea.Time < mea.Time(end) - mea.Padding(2));
		end
		function center = get.Center(mea)
			[~, center] = min(sum((mea.Position - mean(mea.Position)).^2, 2));
		end
		function pr = get.PeakRatio(mea)
			pr = mea.PeakRatio;
			if ~isempty(pr), return; end
			pr = numel(mea.get_wave_times('lfp')) / numel(mea.get_wave_times('lfp_lo'));
			mea.PeakRatio = pr;
		end
		function out = compute_cohgram(mea, varargin)
			% out = mea.compute_cohgram(::'morlet'::, ::'delays'::)
			MW = false;
			NODELAYS = true;
			chars = cellfun(@ischar, varargin);
			if ismember('morlet', lower(varargin(chars))), MW = true; end
			if ismember('delays', lower(varargin(chars))), NODELAYS = false; end
			
			T = 1;
			movingwin = [T .01];
			thresh = 5e-2; 
			MIN_RATIO_FINITE = 0.2;
	
			params.err = [1 thresh];
			params.fpass = [1 50];
			params.tapers = [3 T 1];
			params.Fs = mea.SamplingRate;
			data = mea.Data;
			ctr = mea.Coherent;
			
			N = size(data, 2);
			data1 = data(:, ctr) .* ones(1, N-1);
			data2 = data(:, ~(1:N == ctr));
			if MW
				t = mea.Time;
% 				[C,phi, ~, ~, ~, f, ~, ~, confC, ~] = ...
% 					morlet_coherence(data1, data2, params);
				[~, f, ~, W] = mea.morlet_wavelet;
				W = permute(W, [2 1 3]);
				[nt, nf] = size(W, [1 2]);
				phi = angle(W);
				phi_ctr = phi(:, :, ctr);
				phi(:, :, ctr) = [];
				phi = angle(exp(1j*(phi - phi_ctr)));
				S = abs(W);
				for tt = nt:-1:1
					temp = corrcoef(squeeze(S(tt, :, :)));
					C(tt, :) = temp(ctr, :);
				end
				
				C(:, ctr) = [];
				C = repmat(C, 1, 1, nf);
				C = permute(C, [1 3 2]);
				confC = .8;  % arbitrary... think on this
				
			else
				[C,phi,~,~,~,t,f,confC,~]=cohgramc(data1,data2,movingwin,params);
				 t = t - mea.Padding(1);
			end
			if NODELAYS
				out.C = C; 
				out.phi = phi;
				out.t = t;
				out.f = f;
				out.confC = confC;
				return
			end
			[nt, nf, np] = size(C);
% 			df = mean(diff(f)); 
			winsz = max(3, 3/movingwin(1));  % Hz (this makes sure you get at least 3 points at frequency res)
			transform = @(A) reshape(permute(A, [2 1 3]), nf, []);
			Cf = transform(C);  
			phif = transform(phi); 
			clear C phi
			
			for type = ["group" "phase"]
				delays = dblr(phif, f, type, Cf, confC(1), winsz);
				delaysR = reshape(delays, nf, nt, []);
				clear delays 

				pos = mea.Position(~(1:N == ctr), :);
				warning('off', 'stats:statrobustfit:IterationLimit');
				H = [0 1 0; 0 0 1];  
				c = [0; 0];

				[Z, pdel, pct, V] = deal(nan(nf, nt));
				for ii = 1:nf  % For each frequency
					if ~mod(ii, 100), fprintf('ii=%d/%d\n', ii, nf), end
					for jj = 1:nt  % ... and time point
						delays2fit = squeeze(delaysR(ii, jj, :));  % ... collect the delays for each pair
						num_unique = numel(unique(delays2fit(isfinite(delays2fit))));  % ... calculate the number of finite unique values
						if num_unique < 3, continue; end  % check that there are enough unique points to fit a plane
						finite = numfinite(delays2fit);  % Calculate the number of finite delays
						pct(ii, jj) = finite / np;  % (FYI only)
						if finite <= max(MIN_RATIO_FINITE * np, 3); continue; end  % check enough delay data is not NaN.
						[beta,stats] = robustfit(pos, delays2fit, 'fair');  % fit the delay vs two-dimensional positions
						ptemp = linhyptest(beta, stats.covb, c, H, stats.dfe);  % Compute significance of fit
						if ptemp >= thresh || isnan(ptemp); continue; end  % Stop if fit is not significant
						Vt = pinv(beta(2:3));  % velocity is the pseudoinvervse of the fitted wave
						Z(ii, jj) = angle([1 1i] * Vt(:));  % Z is the angle of the velocity
						V(ii, jj) = abs([1 1i] * Vt(:));
						pdel(ii, jj) = ptemp;  % Store p-value
					end
				end
				
				out.(['Z' char(type)]) = Z;
				out.(['V' char(type)]) = V;
				out.(['pct_' char(type)]) = pct;
				out.(['p_' char(type)]) = pdel;
			end
		
			out.t = t;
			out.f = f;
			
		end
		
		function S = get.Spectrum(mea)
			S = mea.Spectrum;
			if isempty(S)
				[ss, frq] = mea.morlet_wavelet;
				S.power = ss';
				S.frq = frq;
				mea.Spectrum = S;
			end
		end
		function out = exclude_by_coh(mea, electrode)
			% Excludes channels that have coh<.8 relative to electrode
			% Usage:
			%     out = mea.exclude_by_coh(electrode=mea.Coherent)
			if nargin < 2, electrode = mea.Coherent; end
			t_inds = (mea.Time > 0) & (mea.Time < mea.Time(end) - mea.Padding(2));
			data = mea.Data(t_inds, :);
			out.electrode = electrode;
			[cxy, f] = mscohere(data(:, electrode), ...
				data(:, [1:electrode-1, electrode+1:end]), ...
				[], [], [], mea.SamplingRate);
			out.cxy = cxy;
			out.f = f;
			mask = f < 13;
			ch = 1:size(data, 2);
			exclude = find(mean(cxy(mask, :)) < .8);
			exclude(exclude >= electrode) = exclude(exclude >= electrode) + 1;
			ch(exclude) = [];
			mea.Coherent = find(ch == electrode);
			mea.ExcludedChannels = exclude;
		end
		
		function out = get.PairwiseCoherence(mea)
			out = mea.PairwiseCoherence;
			if ~isempty(out), return; end
			t_inds = (mea.Time > 0) & (mea.Time < mea.Time(end) - mea.Padding(2));
			data = mea.Data(t_inds, :);
			N = size(data, 2);
			pairs = nchoosek(1:N, 2);
			[cxy, f] = ...
				mscohere(data(:, pairs(:, 1)), data(:, pairs(:, 2)), ...
				[], [], [], mea.SamplingRate);
			out.cxy = cxy;
			out.f = f;
			out.pairs = pairs;
		end
		function ctr = get.Coherent(mea)
			ctr = mea.Coherent;
			if ~isempty(ctr), return; end
			out = mea.PairwiseCoherence;
			mask = out.f < 13; 
			value = mean(out.cxy(mask, :));
			mat = sparse(out.pairs(:, 1), out.pairs(:, 2), double(value));
			is_coh = full(mat) > .8;
			num_paired = sum(is_coh, 1) + [0; sum(is_coh, 2)]';
			[~, ctr] = max(num_paired);
			mea.Coherent = ctr;
		end
		
		function out = cluster_by_coh(mea)
			
			cxy = mea.PairwiseCoherence;
			[out.coeff, out.score, out.latent, ~, out.explained] = ...
				pca(cxy.cxy', 'NumComponents', 5);
			
		end
		function raw = get.Raw(mea)
			raw = mea.Raw;
			if isempty(raw); mea.load; raw = mea.Raw; end
			raw(:, mea.BadChannels) = [];
			raw(:, mea.ExcludedChannels) = [];
		end
% 		function set.Raw(mea, raw)
% 			mea.Raw(:, mea.locs) = raw;
% 		end
		function p = get.Path(mea)
			file = dir(mea.Path);
			p = fullfile(file(1).folder, file(1).name);
		end
		function set.Path(mea, value)
			file = dir(value);
			mea.Path = fullfile(file.folder, file.name);
		end
		
		function mov = preview(mea, times, ax, varargin)
			% mov = preview(times, ax=axes(figure), ::imagesc directives:: )
			if nargin < 3 || isempty(ax), ax = axes(figure); end
			t_inds = mea.time2inds(times);
			N = length(t_inds);
			im = imagesc(ax, nan(max(mea.Position)), varargin{:});
			ttl = title(ax, '0');
			mov(N) = getframe(ax);
			for ii = 1:N
				im.CData(mea.locs) = mea.Data(t_inds(ii), :);
				set(ttl, 'string', num2str(times(ii)));
				drawnow;
				mov(ii) = getframe(ax);
			end
		end
		function time = get.Time(mea)
			time = downsample(mea.AllTime(), mea.skipfactor);
		end
		function time = get.AllTime(mea)
			time = mea.AllTime();
		end
		
		function set.SamplingRate(mea, value)
			mea.SamplingRate = value;
		end
		function SR = get.SamplingRate(mea)
			SR = min(max(mea.SamplingRate, 1), mea.SRO);
		end
		function set.skipfactor(mea, value)
			if isempty(value), return, end
			mea.skipfactor = max(min(value, mea.SRO), 1);
			mea.SamplingRate = mea.SRO / mea.skipfactor;
		end
		
		function skipfactor = get.skipfactor(mea)
			skipfactor = round(mea.SRO / mea.SamplingRate);
		end

		function data = get.Data(mea)
% 			data = single(downsample(mea.Raw, mea.skipfactor));
			data = single(resample(double(mea.Raw), 1, mea.skipfactor));
		end
		
		function lfp = get.lfp(mea)
			lfp = mea.lfp;
			if isempty(lfp)
				fprintf('Filtering lfp ... ')
				lfp = mea.filter(mea.Data, mea.SamplingRate, [1 50]); 
				mea.lfp = lfp;
				fprintf('Done. \n')
			end
		end

		function artefacts = get.artefacts(mea)
			if isempty(mea.artefacts), mea.mua; end
			artefacts = mea.artefacts;
		end
		function mua = get.mua(mea)
			mua = mea.mua;
			if isempty(mua)
				disp('Filtering mua ...')
				mua = mea.filter(mea.Raw, mea.SRO, [300 3000]);
				
				disp('Done')
				mea.artefacts = mea.get_artefacts(mua, mea.params.artefact_thresh);
				mua(mea.artefacts) = nan;
				mea.mua = mua;
			end
		end
		function EI = get.event_inds(mea)
			EI = mea.event_inds;
			
			if isempty(EI)
				data = mea.baseline_zscore(mea.mua, mea.AllTime);
				intervalM = mea.SRO / 1e3 * mea.params.min_dist;  % samples per MIN_DIST

				% find events on each channel
				events = false(size(data));
				for ch = 1:size(data, 2)
					temp = mea.interp_nan(data(:, ch));
					if ~any(abs(zscore(temp)) > mea.params.event_thresh)
						continue
					end
					[~, inds] = findpeaks(-temp, ...
						'minpeakdistance', intervalM, ...
						'minpeakheight', mea.params.event_thresh);  
					events(inds, ch) = true;
				end
				EI = sparse(events);
				mea.event_inds = EI;
			elseif ~issparse(EI)
				sz = size(mea.Raw);
				[i, j] = ind2sub(sz, EI);
				EI = sparse(i, j, true, sz(1), sz(2));
			end
			
% 			if ~isempty(EI) && ~iscell(EI)
% 				EI = {EI, size(mea.Raw)};
% 			end
		end
		function ET = get.event_times(mea)
			ET = mea.event_times;
			if isempty(ET)
				EI = mea.event_inds;
				[r, ch] = find(EI);
				times = mea.AllTime(r);
				ET = {ch(:), times(:)};
				mea.event_times = ET;
			end
		end
		
		function events = get.mua_events(mea)
			events = mea.event_inds;
		end
		
		function lfp_lo = get.lfp_lo(mea)
			lfp_lo = mea.lfp_lo;
			if isempty(lfp_lo)
				fprintf('Filtering lfp_lo ... ')
				lfp_lo = mea.filter(mea.Data, mea.SamplingRate, [1 13]);
				fprintf('Done. \n')
				mea.lfp_lo = lfp_lo;
			end
		end
		
		function pos = get.Position(mea)
			pos = mea.Position;
			pos(mea.BadChannels, :) = [];
			pos(mea.ExcludedChannels, :) = [];
		end
		
			
		function fr = get.firing_rate(mea)
			if isempty(mea.firing_rate)
				fr = mea.compute_firing_rate(mea); 
				mea.firing_rate = fr;
			else
				fr = mea.firing_rate;
			end
		end
		function m = all_metrics(s)
			m = fieldnames(s.Fits);
		end
		function m = get.metrics(s)
			m = s.metrics;
			if isempty(m)
				m = s.all_metrics;
				s.metrics = m;
			end
		end
		
		function bsi = get.BSI(s)
			bsi = s.BSI;
			if isempty(s.BSI)
				bsi = s.compute_bsi(s);
				s.BSI = bsi;
			end
			
		end
		function wt = get.WaveTimes(s)
			if isempty(s.WaveTimes)
				wt = s.compute_wave_times(s);
				s.WaveTimes = wt;
			else
				wt = s.WaveTimes;
			end
		end
		function [inds, dist] = time2inds(mea, times)
			inds = interp1(mea.Time, 1:length(mea.Time), times, 'nearest', 'extrap');
			dist = times - mea.Time(inds);
		end
		function locs = get.locs(s)
			locs = sub2ind(max(s.Position), s.Position(:, 1), s.Position(:, 2));
		end
		function exclude_by_std(s)
			sd = isoutlier(std(s.Data));
			ch = 1:100;  % assuming all meas have <=100 channels
			ch(s.BadChannels) = [];
			s.BadChannels = unique([ch(sd) s.BadChannels]);
		end
		function excluded = exclude_channels(mea, excluded)
			EXCLUDE = mea.params.exclude_rate;  % (default: 6 spikes/sec)
			REQUIRE_INCREASE = mea.params.require_increase;  % (default: true)

			seizure_inds = mea.Time > 0 & mea.Time < (mea.Time(end) - mea.Padding(2)); % Indicate when seizure is occurring
			if nargin < 2
				excluded = find(mean(mea.firing_rate(seizure_inds, :)) < EXCLUDE);            % Only use channels with mean firing rate at least 6 spikes per second

				if REQUIRE_INCREASE
					% Find channels that don't have significant increases in firing rate during 
					fr = mea.firing_rate;
					mn = mean(fr(mea.Time < 0, :));    % get mean firing rate at baseline
					sd = std(fr);                  % ... and standard deviation over full recording
					fr = (fr - mn) ./ sd;          % normalize smoothed firing rate

					% the first half of the seizure
					halfWay = (mea.Time(end) - mea.Padding(2)) / 2;
					excluded = ...
						unique([ ...
							excluded ...
							find(max(fr(mea.Time < halfWay, :)) < 2, 1) ...
						])';

				end
			end
			ch = 1:max([mea.BadChannels(:); excluded(:)]);
			ch(mea.BadChannels) = [];
			mea.BadChannels = unique([mea.BadChannels(:); ch(excluded)']);
			
		end
		
		function inds = which_t(s, t0)
			[~, inds] = min(abs(s.Time - t0));
		end
		
	end
	
	methods % wave fitting
		function F = flow(s)
			F = WaveFlow(s);
		end

		function M = max_descent(mea, times, varargin)
			N = numel(times);
			M(N) = MaxDescent(varargin{:});
% 			mea.MaxDescentData = zscore(mea.filter(mea.Data, mea.SamplingRate, M(N).FBand)); 
			for ii = 1:N
				t0 = times(ii);
				M(ii) = MaxDescent(mea, t0, varargin{:});		
			end
			M = WaveProp.resize_obj(M);
			M.Name = mea.Name;
		end
		
		function E = dir_events(mea, times, varargin)
			N = numel(times);
			E(N) = Events;
			for ii = 1:N
				t0 = times(ii);
				E(ii) = Events(mea, t0, varargin{:});		
			end
			E = WaveProp.resize_obj(E);
			E.Name = mea.Name;
		end
		function D = delays(mea, times, varargin)
			N = numel(times);
			wng = warning;
			warning('off');
			D(N) = Delays(mea, times(N), varargin{:});
			num_workers = getenv('NSLOTS');
			if 1  %isempty(num_workers)
				
				tic
				for ii = 1:N-1
					if ~mod(ii, 10)
						fprintf('Computed delay %d/%d (%0.2f seconds)\n', ii-1, N, toc); 
						tic; 
					end
					t0 = times(ii);
					D(ii) = Delays(mea, t0, varargin{:});		
				end
			else
% 				clust = parcluster;
% 				parpool(clust);				
				delay_fun = @(t0) Delays(mea, t0, varargin{:});
				parfor(ii = 1:N, str2double(num_workers))
					warning('off');
					t0 = times(ii);
% 					D(ii) = Delays(mea, t0, varargin{:});
					D(ii) = delay_fun(t0);
					if ~mod(ii, 10)
						fprintf('Computed delay %d/%d (parfor)\n', ii-1, N); 
					end
				end
% 				delete(gcp('nocreate'));
			end
			D = WaveProp.resize_obj(D);
			D.Name = mea.Name;
			warning(wng);
			
		end
	
		function [wt, peaks] = get_wave_times(s, method, min_peak_height)
			% [wt, peaks] = get_wave_times(method='events', min_peak_height=-Inf)
			if nargin < 2, method = 'events'; end
			if nargin < 3, min_peak_height = -Inf; end
			switch method
				case 'bsi'
					data = zscore(mean(s.BSI, 2));
					minProm = 0;
				case 'lfp'
					data = -zscore(mean(s.lfp, 2));
					minProm = 0;
				case 'lfp_lo'
					data = -zscore(mean(s.lfp_lo, 2));
					minProm = 0;
				case 'events'
					data = mean(s.firing_rate, 2);
					data = data ./ nanstd(data);
					minProm = 0;
				otherwise
					error('Method not recognized')
			end
			[pks, lcs, w, p] = findpeaks(data, s.SamplingRate, ...
				'MinPeakHeight', min_peak_height, ...
				'MinPeakDistance', 10*1e-3, ...
				'MinPeakProminence', minProm);
% 				'MinPeakHeight', 1);
			lcs = round(lcs * s.SamplingRate);
			outliers = isoutlier(pks);
% 			wt = s.Time(locs(~outliers));
			wt = s.Time(lcs);
			peaks = struct('pks', pks, 'locs', lcs, ...
				'w', w, 'p', p, 'outliers', outliers);
			
			
		end
		
		function [S, fs] = pspectrum(mea, data, fs, range)
			% [S, fs] = pspectrum(data=Data, fs=SamplingRate, range=[0 100])
			if nargin < 4, range = [0 100]; end
			if nargin < 3, fs = mea.SamplingRate; end
			if nargin < 2, data = mea.Data; end
			
			[S, fs] = pspectrum(data, fs);
			S = S(fs > range(1) & fs < range(2), :);
			fs = fs(fs > range(1) & fs < range(2));
			S = pow2db(S.*fs.^2);
		end
		function noisy_sig = add_noise(mea, snr_dB)
			signal = mea.Data;
			brown_noise = randnd(-2, [length(signal), max(mea.Position)]);  % brownian noise: -2; pink noise: -1
			noise = brown_noise(:, mea.locs);
			signal_power = rms(signal);
			noise_power = rms(noise);
			scale_factor = signal_power ./ (10^(log10(snr_dB)/2) .* noise_power);
			noisy_sig = signal + scale_factor .* noise;
			
		end
		function out = change_points(mea, data, thresh, min_dist)
			% out = mea.change_points(data='fr',
			% thresh=30*mea.SamplingRate, min_dist=mea.SamplingRate)
			if nargin < 4 || isempty(min_dist), min_dist = 1; end
			if nargin < 3 || isempty(thresh), thresh = 30; end
			if nargin < 2 || isempty(data) 
				data = 'fr';
			end
			if ischar(data)
				data = validatestring(data, ...
					{'fr', 'firing', 'spectrum', 'coherence', 'direction', ...
					'mean_fr' 'isi', 'peak_height', 'lfp'});
				switch data
					case {'fr', 'firing'}
						data = mea.firing_rate; 
						data = zscore(data);
						out.t = mea.Time;
						thresh = thresh*mea.SamplingRate;
						min_dist = min_dist * mea.SamplingRate;
						stat = 'linear';
					case 'spectrum'
						S = mea.Spectrum;
						data = zscore(pow2db(S.power));
						out.frq = S.frq;
						out.t = mea.Time;
						thresh = thresh*mea.SamplingRate;
						min_dist = min_dist * mea.SamplingRate;
						stat = 'mean';
					case 'coherence'
						coh = mea.compute_cohgram;
						coh.C(coh.C < coh.confC(1)) = nan;
						data = squeeze(nanmean(coh.C, 3));
						data(~isfinite(data)) = 0;
						out.t = coh.t;
						out.frq = coh.f;
						data = zscore(data);
						SR = 1 / diff(coh.t(1:2));
						thresh = thresh * SR;
						min_dist = min_dist * SR;
						stat = 'mean';
					case 'direction'
						F = WaveProp.load(mea, 'metrics', {'D1xwh'});
						F = F.D1xwh;
						out.t = F.t0;
						data = [0; diff_phase(F.Direction) ./ diff(F.t0)];
						mask = isfinite(data);
						data = data(mask);
						out.t = out.t(mask);
% 						thresh = nanstd(data);
						thresh = 150;
						min_dist = 2 / median(diff(F.t0));
						stat = 'std';
					case 'mean_fr'
						data = mean(mea.firing_rate, 2);
						out.t = mea.Time;
						thresh = thresh/2 * mea.SamplingRate;
						min_dist = min_dist * mea.SamplingRate;
						stat = 'rms';
					case {'isi', 'peak_height'}
						lfp = -zscore(mea.lfp);
						N = size(lfp, 2);
						clear isi pks
						for pp = N:-1:1
							[pks{pp}, tpts{pp}] = ...
								findpeaks(lfp(:, pp), mea.Time, ...
								'MinPeakProminence', 2);
							isi{pp} = diff(tpts{pp});
							t{pp} = tpts{pp}(2:end);
						end
						if strcmpi(data, 'isi')
							[t, so] = sort(cat(2, t{:}), 'ascend');
							data = cat(2, isi{:});
							data = data(so);
							out.t = t;
						else
							[t, so] = sort(cat(2, tpts{:}), 'ascend');
							data = cat(1, pks{:});
							data = data(so);
							out.t = t;
						end
						stat = 'linear';
						min_dist = 10;
					case 'lfp'
						data = -zscore(mea.lfp);
						out.t = mea.Time;
						stat = 'linear';
						thresh = mea.SamplingRate;
						min_dist = min_dist * mea.SamplingRate;
				end
			end
			out.data = data;
			[q, r] = findchangepts(data', ...
				'Statistic', stat, ...
				'MinThreshold', thresh, ...
				'MinDistance', round(min_dist));
% 			[q, r] = findchangepts(data', ...
% 				'Statistic', stat, ...
% 				'MaxNumChanges', 5, ...
% 				'MinDistance', min_dist);
			out.CP1 = q(:)';
			out.r1 = r;
			mask = false(size(data, 1), 1);
			mask(q) = true;
			groups = cumsum(mask) + 1;
			out.group = uint16(groups);
			[q, r] = findchangepts(data', ...
				'Statistic', 'mean', ...
				'MinThreshold', thresh/10, ...
				'MinDistance', round(min_dist));
			out.CP2 = q(:)';
			out.r2 = r;
		end
		function [S, frq, coi, W] = morlet_wavelet(mea, data, fs, range, varargin)
			% [S, frq, coi, W] = mea.morlet_wavelet(data=Data, fs=SamplingRate, range=[1 50], ::'mean'::)
			if nargin < 4 || isempty(range), range = [1 50]; end
			if nargin < 3 || isempty(fs), fs = mea.SamplingRate; end
			if nargin < 2 || isempty(data), data = mea.Data; end
			if ismember('mean', lower(varargin)), data = mean(data, 2); end
			
			if nargout > 3
				[S, frq, coi, W] = ...
					compute_wavelet(data, range(1), range(2), fs);
			else
				[S, frq, coi] = ...
					compute_wavelet(data, range(1), range(2), fs);
			end
		end
		function [outliers, idx, dist, C] = cluster_traces(mea, data, kN, varargin)
			% [outliers, idx, dist, C] = cluster_traces(data=Data, kN=Nch/4, ::'tsne'::)
			if nargin < 3 || isempty(kN), kN = floor(size(mea.Data, 2) / 4); end
			if nargin < 2 || isempty(data), data = mea.pspectrum'; end
			if ismember('tsne', lower(varargin))
				data = tsne(data);
			end
			[idx, C] = kmeans(data, kN);
			dist = sum(squareform(pdist(C)));
			outliers = find(isoutlier(dist));
		end
	end
	
	
	
	methods (Static)
		
		function bsi = compute_bsi(s)
			N = numel(s.event_times);
			edges = [s.Time(:); Inf];
			bsi = nan(size(s.Data));
			for ii = 1:N
				bsi(:, ii) = histcounts(s.event_times{ii}, edges);
			end
		end
		
		function data = baseline_zscore(data, time)
			if ~any(time < 0) 
				mn = 0; 
				sd = .01 * range(data);
			else
				mn = mean(data(time < 0, :), 1, 'omitnan');     % get the mean of the preictal baseline
				sd = std(data(time < 0, :), 'omitnan');            % ... and the sd
			end
			data = (data - mn) ./ sd;                       % zscore based on preictal state

		end
		function artefacts = get_artefacts(data, thresh)
			artefacts = (abs(zscore(data)) > thresh); 
		end
		function data = interp_nan(data)
			if ~any(isnan(data)), return, end
			x = 1:length(data);
			mask = isfinite(data);
			data = interp1(x(mask), data(mask), x, 'nearest', 'extrap');
		end
		
		function [patient, seizure] = get_info(path)
			[~, data] = fileparts(path);
			data = strsplit(data, '_');
			patient = data{1};
			seizure = data{2}(8:end);
		end
		
		
		function data = filter(data, sampling_rate, band)
        % data = filter(data, sampling_rate, band);  % (static)
			if band(2) > sampling_rate / 2 / 1.1
				band(2) = (sampling_rate / 2 - 1) / 1.1;
				fprintf('Setting CutoffFrequency2 to %.1f\n', band(2)); 
			end
% 			bpFilt = designfilt('bandpassiir','FilterOrder',150, ...
% 				'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2), ...
% 				'SampleRate', sampling_rate);
			if band(1) == 0
				bpFilt = designfilt('lowpassfir', 'filterorder', 150, ...
					'passbandfrequency', band(2), ...
					'stopbandfrequency', 1.1*band(2), ...
					'samplerate', sampling_rate);
% 				bpFilt = designfilt('lowpassiir', ...
% 					'PassbandFrequency', band(2), ...
% 					'StopbandFrequency', 1.1*band(2), ...
% 					'PassbandRipple', .5, 'StopbandAttenuation', 65, ...
% 					'SampleRate', sampling_rate);
			else
				bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
					'CutoffFrequency1',band(1),'Cutofffrequency2',band(2), ...
					'SampleRate', sampling_rate);
% 				bpFilt = designfilt('bandpassiir', ...
% 					'StopbandFrequency1', .1 * band(1), 'PassbandFrequency1', band(1), ...
% 					'PassbandFrequency2', band(2), 'StopbandFrequency2', 1.1*band(2), ...
% 					'StopbandAttenuation1', 65, 'PassbandRipple', .5, ...
% 					'StopbandAttenuation2', 65, 'SampleRate', sampling_rate, ...
% 					'MatchExactly', 'passband');
			end
			data = single(filtfilt(bpFilt, double(data)));
		end
		function firing_rate = compute_firing_rate(mea)
			
			samplingRateMs = mea.SRO / 1e3;  % samples per ms
			window = mea.params.fr_window * samplingRateMs;  % number of samples to use in the window

			% Compute spike rate
			firing_rate = downsample(...
					smoothdata(mea.mua_events, 1, 'movmean', window), ...
					mea.skipfactor ...
				) * mea.SRO;  % convert to spikes per second 
		end
		
	end
	
	methods  % Plotting
		function plot_panels(mea, times, layout, h)
			if nargin < 4, h = gcf; end
			if nargin < 3, layout = []; end
			if isempty(layout)
				layout = [ceil(numel(times)/7) 7];
			end
			set(0, 'currentfigure', h);
			inds = mea.time2inds(times);
			data = zscore(mea.Data);
			r = layout(1); c = layout(2);
			temp = nan(max(mea.Position));
			for ii = 1:length(inds)
				subplot(r, c, ii)
				temp(mea.locs) = data(inds(ii), :);
				imagesc(temp, [-2 2]);
				xticks([]);
				yticks([]);
				colormap gray
				if ii == 1
					title(num2str(times(ii)));
				end
			end
			ttl = title(num2str(times(ii)));
			lbl = xlabel(num2str(times(ii)));
			ttl.Position(2:3) = lbl.Position(2:3);
			delete(lbl);
			ttl.VerticalAlignment = 'top';
		end
		function ax = plot_std(mea, window, ax)
			% ax = plot_std(window=1, ax=gca)
			if nargin < 3, ax = gca; end
			if nargin < 2 || isempty(window), window = 1; end
			for M = string(mea.metrics(:)')
				F.(M) = mea.Fits.(M);
			end
			[Z, t0] = WaveProp.WP2mat(F, 'Direction');
			Zstd = circ_std(Z');
			Zstd = smoothdata(Zstd', 1, 'gaussian', window, 'SamplePoints', t0);
			ZstdZ = (Zstd - nanmean(Zstd)) / nanstd(Zstd);
			imagesc(ax, t0, 1, ZstdZ', [-2 2])
			xticklabels(ax, []);
			yticklabels(ax, []);
			yticks(ax, []);
			ax.Tag = ['figs/' mea.patient '_' mea.seizure '_std'];
		end
		function fr = plot(s, ax)
			if nargin < 2, ax = gca; end
			cmap = ax.ColorOrder;
			set(ax, 'nextplot', 'replaceChildren', 'xlim', [s.Time(1), s.Time(end)]);
			frline = gobjects(3, 1);
			dir = gobjects(numel(s.metrics), 1);
			NP = ax.NextPlot;
			
			frline(1) = plot(ax, s.Time, zscore(mean(s.lfp_lo, 2)), ...
				'color', (cmap(3, :) + .8)/2, ...
				'linewidth', 2, ...
				'displayname', 'lfp_lo');
			ax.NextPlot = 'add';
			frline(2) = plot(ax, s.Time, zscore(mean(s.lfp, 2)), ...
				'color', (cmap(1, :) + .8)/2, ...
				'linewidth', 2, ...
				'displayname', 'lfp');
			mn_fr = nanmean(s.firing_rate, 2);
			frline(3) = plot(ax, s.Time, mn_fr ./ nanstd(mn_fr) - pi, ...
				'color', (cmap(2, :) + .8)/2, ...
				'linewidth', 2, ...
				'displayname', 'FR');
			set(ax, 'colororderindex', 1);
			for ii = 1:numel(s.metrics)
				m = s.metrics{ii};
				if strcmp(m, 'F'), sz = 3; else, sz = 8; end
				dir(ii) = scatter(ax, s.Fits.(m).time, s.Fits.(m).Direction, ...
					sz, 'filled', 'MarkerFaceAlpha', .8, ...
					'DisplayName', m);
			end
			ylim(ax, [-pi pi]);
			yticks(ax, -pi:pi/4:pi)
			yticklabels(ax, {'-\pi', '', '-\pi/2', '', '0', ...
				'', '\pi/2', '', '\pi'});
			ylabel(ax, {'Direction', 'Normalized Signal'})
			xlabel(ax, 'Time (s)');
			title(ax, [s.patient s.seizure]);
			ax.NextPlot = NP;
			fr.ax = ax;
			fr.frline = frline;
			fr.dir = dir;
			fr.lgd = legend(ax, dir);
			fr.name = ['figs/' s.patient '_' s.seizure '_plot'];
		end
		function imagesc(s, data, t0, ax)
			if nargin < 4, ax = gca; end
			if nargin < 3 || isempty(t0), t0 = 0; end
			if nargin < 2 || isempty(data), data = mea.Data; end
			temp = nan(max(s.Position));
			ind = s.time2inds(t0);
			temp(s.locs) = data(ind, :);
			imagesc(ax, temp);
		end
	end
	
end



