classdef IW < handle
    
    properties
        
        name
        fmax = 2
        wave = 1
        MinPeakHeight = 2  % Defines the threshold for candidate IW crossing times in standard deviations
        MinPeakFr = 30
        Method char {mustBeMember(Method, {'fr', 'pwr', 'frg'})} = 'fr'  % 'pwr' or 'fr'
        DiffsOrPeaks = 'peaks'
        FiringRateWin = 20  % in ms
        SmoothingArgs = {'movmean', 1}  % alternately, use {'gaussian', .2}; {method, seconds}
        
        MaxTemplates = inf
        W  % Big smoothing window for looking for IW
        MinElectrodes = 10
    end
    
    properties (Transient = true)
        ManualOuts
        all_pks   % Used in get_stats plots
        all_locs
        all_durs

        mea
        iw_templates
        locs
        fr
        fr_movmean
        fr_gaussian
        fr_smooth
        fr_at_peak
        iw_fwhm
    end
    
    properties (Dependent = true)
        time
        nch
        pks
        t_mask
        V
        speed
        phi
        num_waves
        center
        range
        show
        position
        
        onsets
        durs
        outliers
        GridSize
    end
    
    methods
        
        function self = IW(mea, varargin)
            for ii = 1:2:nargin-1
                field = validatestring(varargin{ii}, properties(self));
                self.(field) = varargin{ii + 1};
            end
            if isnumeric(mea), mea = MEA(mea); end
            self.mea = mea;
            self.name = mea.Name;
            if ~isempty(self.FiringRateWin)
                self.mea.params.fr_window = self.FiringRateWin;
            end

        end
        
        
        % basic mea properties
        function nch = get.nch(self); nch = size(self.mea.Data, 2); end
        
        function locs = get.locs(self) 
            if isempty(self.locs), self.locs = self.mea.locs; end
            locs = self.locs;
        end
        
        function gs = get.GridSize(self)
            gs = self.mea.GridSize;
        end
        
        function set.GridSize(self, value)
            self.mea.GridSize = value;
        end
                
        function mea = get.mea(self)
            if isempty(self.mea)
                pat = strsplit(self.name, '_');
                self.mea = MEA(sprintf('%s/%s_Neuroport_10_10.mat', ...
                    pat{1}, self.name));
                self.mea.params.fr_window = self.FiringRateWin;
            end
            mea = self.mea;
        end
        
        function fr = get.fr(self)
            if isempty(self.fr) || self.FiringRateWin ~= self.mea.params.fr_window
                self.mea.firing_rate = [];
                self.mea.params.fr_window = self.FiringRateWin;
                self.fr = self.mea.firing_rate;
                self.fr_smooth = [];
            end
            fr = self.fr;
        end
        
        function t = get.time(self); t = self.mea.Time; end
        
        
        
        
        %% Get IW templates
        
        
        
        function tpl = get.iw_templates(self)
            
            RECOMPUTE = isempty(self.iw_templates) ...
                || ~strcmpi(self.iw_templates.method, self.Method) ...
                || ~strcmpi(self.DiffsOrPeaks, self.iw_templates.mdorpeaks);
            
            if RECOMPUTE
                self.compute_IW_templates;
            end
            tpl = self.iw_templates;
        end
        
        function ind = main_wave(self)
            % Choose the wave with the highest firing. Set
            % non-participating electrodes to 0 and use the mean to compute
            % the firing.
            fr_ = mean(fillmissing(self.iw_templates.firing_rate, 'constant', 0), [2 3]);
            [~, ind] = max(fr_);
        end
        
        function M = max_descent_IW0(iw, times, varargin)
            mea_ = iw.mea;
            MIN_FR = iw.MinPeakFr;
            
            
            data = max(iw.fr_smooth, iw.MinPeakFr) - iw.MinPeakFr;
%             data = iw.fr_smooth;
            dataN = normalize(data, 'scale');
            
            S = warning; warning off;
			N = numel(times);
			M(N) = MaxDescent(varargin{:});
            D0 = mea_.MaxDescentData;
            dataN(data < MIN_FR) = 0;  % Require min firing rate of MEN_FR
            mea_.MaxDescentData = -dataN;
%             M = MaxDescent(mea_, times, ...
%                 'halfwin', 1, 'diffsorpeaks', iw.DiffsOrPeaks, ...  % defaults
%                 varargin{:});
			for ii = 1:N
				t0 = times(ii);
				M(ii) = MaxDescent(mea_, t0, ...
                    'halfwin', 1, ...
                    'diffsorpeaks', iw.DiffsOrPeaks, ...
                    varargin{:});		
			end
			M = WaveProp.resize_obj(M);
            M.MinFinite = 10;  % Set to 10 to accommodate slow movement relative to window
			M.Name = ['IW_' mea_.Name];
            
            warning(S);
            mea_.MaxDescentData = D0;
        end
        
        function [M, fwhm] = max_descent_IW(iw, times, varargin)
            % Use max descent method to compute IW-like wave propagation
            
            mea_ = iw.mea;
            time_ = iw.time;
            MIN_FR = iw.MinPeakFr;
            
            data = iw.fr_smooth;
            dataN = normalize(max(data, MIN_FR), 'scale');
            
            S = warning; warning off;
            D0 = mea_.MaxDescentData;
            mea_.MaxDescentData = -dataN;
            M = MaxDescent(mea_, times, ...
                'halfwin', 1, 'diffsorpeaks', iw.DiffsOrPeaks, ...  % defaults
                varargin{:});

            M.TOA = M.use_largest_cluster(M.TOA);
            M.MinFinite = 10;  % Set to 10 to accommodate slow movement relative to window
			M.Name = ['IW_' mea_.Name];
            
            warning(S);
            mea_.MaxDescentData = D0;
            
            % Get the IW width (duration) of the IW on each electrode
            dat = M.TOA + M.t0 - M.HalfWin;
            fwhm = zeros(size(dat));
            
            [~, locs_w, ww] = arrayfun(@(ii) findpeaks(dataN(:, ii), time_), 1:size(dataN, 2), 'uni', 0);
            for ii = 1:numel(ww)
                li = mea_.locs(ii);
                temp_tq = dat(:, li);
                temp_t_all = locs_w{ii}(:);
                [~, nearest] = min(abs(temp_tq - temp_t_all'), [], 2);
                temp_w = ww{ii}(nearest);
                
                % don't add data where there isn't any, or where the times
                % don't match
                t_diff = temp_tq - temp_t_all(nearest);
                mask = isfinite(t_diff) & abs(t_diff) < 3e-3;  % within 3 ms
                temp_w(~mask) = nan;
                
                fwhm(:, li) = temp_w;
            end

            
            
        end

        function [M, s] = save_IW_fits(iw, method, save_bool)
            % Creates the IW files in WaveFits
            if nargin < 2, method = iw.Method; end
            if nargin < 3, save_bool = true; end
            s = [];
            iw.Method = method;
            
            switch method
                case {'fr', 'frg', 'pwr'}
                    M = iw.max_descent_IW(iw.fr_smooth, ...
                        'diffsorpeaks', 'peaks');
                case 'pwr_md'
                    M = iw.max_descent_IW(iw.fr_smooth, ...
                        'diffsorpeaks', 'diffs');
            end
            Mname = ['Miw_' method];
            if save_bool
                s = mkdir(['WaveFits/' iw.mea.Name]);
                save(['WaveFits/' iw.mea.Name filesep Mname], 'M');
            end
        end
        
        
        function out = compute_IW_templates(iw, method, varargin)
            % [out, M] = compute_IW_templates(iw, method, varargin)
            
            max_templates = iw.MaxTemplates;  % Return at most this many templates
            win = iw.W;  % Use this as the min peak window (otherwise, autodetect)
            mea_ = iw.mea;
            time_ = mea_.Time;
            MIN_FR = iw.MinPeakFr;  % channels must hit this firing rate to be considered in IW

            if nargin < 2 || isempty(method), method = iw.Method; end
            assert(ischar(method), 'Use the <method> so that the right firing rate is used')

            % Get the firing rate according to the indicated method
            iw.Method = method;
            fr_ = iw.fr_smooth;


            % Normalize the firing rate based on std
            frN = normalize(fr_, 'scale');
            
            % Auto-detect min peak window length: Use the width of the
            % largest peak as the window to differentiate between peaks
            if isempty(win)  
                num_hi_fr = movmean(sum(frN > iw.MinPeakHeight, 2), ...
                    mea_.SamplingRate*1);  % Use a 1 s smoothing window since this should wash out fast changes
                [~, ~, win] = findpeaks(num_hi_fr, time_, ...
                    'SortStr', 'descend', 'Npeaks', 1, ...
                    'WidthReference', 'halfprom');  % halfprom is what you want here... don't change this! (halfheight sounds like what you want but it truncates; prom is what you want)
                if isempty(win), win = 1; end
                win = max(min(win, 60), 2);
            end

            

            % Find time of peaks on each channel
            [~, locs_t] = arrayfun( ...
                @(ii) findpeaks(frN(:, ii), time_, ...
                'minpeakheight', iw.MinPeakHeight, 'minpeakdistance', win), ...
                1:size(fr_, 2), 'uni', 0);
            ch = arrayfun(@(ii) ii*ones(size(locs_t{ii})), 1:size(fr_, 2), 'uni', 0);

            
            % Convert peak times to locs
            locs_t = cat(2, locs_t{:});
            ch = cat(2, ch{:});
            [~, locs_i] = min(abs(time_(:) - locs_t));

            
            % Create a matrix of peaks and determine how many channels have
            % a peak in each sliding window
            [m, n] = size(fr_);
            temp = movmax(full(sparse(locs_i, ch, 1, m, n)), win * mea_.SamplingRate);
            N = sum(temp, 2);
            
            
            % Find time points where many electrodes have a peak
            [~, times, durs_, proms_] = findpeaks(N, time_, ...
                'MinPeakDistance', win, 'MinPeakHeight', iw.MinElectrodes);

            
            % Compute max descent (or min peak) on time points
            [M, fwhm] = iw.max_descent_IW(times, 'halfwin', win/2, varargin{:});
            dat = WaveProp.use_largest_cluster(M.Data) + M.time - M.HalfWin;
            [~, inds] = min(abs(time_(:) - times));
            fr_max = movmax(fr_, win * mea_.SamplingRate); 
            
            
            fr_max = fr_max(inds, :);
            fr_max(isnan(dat(:, mea_.locs)) | fr_max < MIN_FR) = nan;
            fr_maxN = normalize(fr_max, 2);
            fr_max(fr_maxN < -2) = nan;  % exclude any channels with low firing rates
            fr_max = mea_.make_3d(fr_max);
            dat(isnan(fr_max)) = nan;
            fr_max(isnan(dat)) = nan;
            fwhm(isnan(dat)) = nan;
            
            
            % Eliminate times where not enough electrodes are active
            mask = sum(isfinite(dat), [2 3]) >= iw.MinElectrodes;

            
            % Return the highest firing templates
            if max_templates < sum(mask)
                fr_max(~mask, :, :) = 0;
                [~, so] = sort(sum(fillmissing(fr_max, 'constant', 0), [2 3]), 'descend');
                mask(so(max_templates+1:end)) = false;
            end
            
            
            out = struct( ...
                'template', dat(mask, :, :), ...
                'firing_rate', fr_max(mask, :, :), ...
                'fwhm', fwhm(mask, :, :), ...
                'time', M.time(mask), ...
                'win', win, ...
                'durs', durs_(mask), ...
                'proms', proms_(mask), ...
                'method', iw.Method, ...
                'mdorpeaks', iw.DiffsOrPeaks, ...
                'M', M);
            
            iw.iw_templates = out;



        end
        

        
        %%
        
        function iw_fwhm_ = get.iw_fwhm(self)
            ind = self.wave;
            tpl = self.iw_templates;
            iw_fwhm_ = tpl.fwhm(ind, self.locs);
        end
        
        function fr = get.fr_at_peak(self)
            ind = self.wave;
            tpl = self.iw_templates;
            fr = tpl.firing_rate(ind, self.locs);
        end
        function onsets = get.onsets(self)
            onsets = self.iw_templates.template(self.wave, self.locs);
        end
        function fwhm = fr_fwhm(self)
            fwhm = self.iw_templates.durs(self.wave);
        end
        function durs = get.durs(self)
            % This used to be the duration of the peak on each channel; now
            % it is the duration of the iw wave as seen on all channels 
            durs = self.iw_templates.durs(self.wave);
        end
        function set.outliers(self, value)
            self.ManualOuts = value;
        end
        function outliers = get.outliers(self)
            outliers = ~isfinite(self.onsets);
        end
           
        
        function V = get.V(self)
            [V, ~] = self.regress;
%             if p >= 0.05, V = [nan nan]; end
        end
        function speed = get.speed(self)
            fit = self.wave_fit;
            speed = fit.speed;
        end
        function phi = get.phi(self)
            phi = atan2(self.V(2), self.V(1));
        end
        
        function num_waves = get.num_waves(self)
            num_waves = numel(self.iw_templates.time);
        end
        function center = get.center(self)
            center = nanmedian(self.onsets(~self.outliers));
        end
        function range = get.range(self)
            % Because the iw starts so early in some cases, I only used the
            % first to last peak times rather than including the halfwidth
            % info
            range = quantile(self.onsets(~self.outliers), [0 1]);
            
%             r0 = self.onsets + [-.5; .5] .* self.durs;
%             r0(:, self.outliers) = [];
%             range = quantile(r0(:), [0 1]);
            
        end
        function pos = get.position(self)
            [p1, p2] = ind2sub(self.mea.GridSize, self.locs);
            pos = [p1, p2];
        end
        
        
        function ap = get.all_pks(self)
            if isempty(self.all_pks)
                [pks_, locs_, durs_] = self.findpeaks([-Inf Inf]);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Visual inspection for multiple waves
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                self.all_pks = cat(1, pks_{:});
                self.all_locs = cat(2, locs_{:});
                self.all_durs = cat(2, durs_{:});
            end
            ap = self.all_pks;
        end
        function al = get.all_locs(self)
            if isempty(self.all_locs), self.all_pks; end
            al = self.all_locs;
        end
        function ad = get.all_durs(self)
            if isempty(self.all_durs), self.all_pks; end
            ad = self.all_durs;
        end

        function frS = get.fr_gaussian(self)
            if isempty(self.fr_gaussian)
                frS = smoothdata(self.fr, 'gaussian', .2, ...
                        'SamplePoints', self.time);
                self.fr_gaussian = frS;
            end
            frS = self.fr_gaussian;
        end
        function frS = get.fr_movmean(self)
            if isempty(self.fr_movmean)
                frS = smoothdata(self.fr, 'movmean', 1, ...
                        'SamplePoints', self.time);
                self.fr_movmean = frS;
            end
            frS = self.fr_movmean;
        end
        function frS = get.fr_smooth(self)
            % Smooth firing rate according to method
            switch self.Method
                case 'frg'
                    frS = self.fr_gaussian;
                case 'fr'
                    frS = self.fr_movmean;
                case 'pwr'
                    frS = self.mea.iw_firing_rate;
                otherwise
                    error('No definition for method <%s>.', self.Method)
            end

        end
        
        
        function [pks, locs, fwhm] = findpeaks(self, t_bounds_)
            if nargin < 2, t_bounds_ = self.t_bounds; end
            S = warning;
            switch self.Method
                case 'pwr'
                    data = self.p_lo;
                    tt = self.t;
                    thresh = self.MinPeakHeight;
                case 'fr'
                    data = self.fr_smooth;
                    tt = self.time;
                    data(tt < 0 | tt > self.seizure_dur, :) = nan;
                    data = normalize(data);
%                     thresh = self.MinPeakFr;
                    thresh = self.MinPeakHeight;
            end
            t_mask_ = tt >= t_bounds_(1) & tt <= t_bounds_(2);
            warning('off')
            [pks, locs, fwhm] = arrayfun( ...
                        @(ich) findpeaks( ...
                        data(t_mask_, ich), tt(t_mask_), ...
                        'minpeakheight', thresh, ...
                        'SortStr', 'descend'), ...
                        1:self.nch, 'uni', 0);
            warning(S)
        end
        function fr_max = peak_fr(self)
            
            t_on = self.onsets ...
                - self.durs/2;
            t_off = self.onsets ...
                + self.durs/2;

            fr_ = self.fr;
            
            [~, ind_on] = min(abs(self.time' - t_on));
            [~, ind_off] = min(abs(self.time' - t_off));
            
            N = size(fr_, 2);
            fr_max = nan(1, N);
            for ich = 1:N
                fr_max(ich) = max(fr_(ind_on(ich):ind_off(ich), ich));
            end

        end
        function fr_mn = mean_fr(self)
            
            t_on = self.onsets ...
                - self.durs/2;
            t_off = self.onsets ...
                + self.durs/2;

            fr_ = self.fr;
            
            [~, ind_on] = min(abs(self.time' - t_on));
            [~, ind_off] = min(abs(self.time' - t_off));
            
            N = size(fr_, 2);
            fr_mn = nan(1, N);
            for ich = 1:N
                fr_mn(ich) = mean(fr_(ind_on(ich):ind_off(ich), ich));
            end

        end
        function sp = prop_sp(self)
            sp = range(self.onsets(~self.outliers)); %#ok<CPROP>
        end
        
        function fit = fit_iw_fwhm(self)
            iw_fwhm_ = self.iw_fwhm;
            fit = fitdist(iw_fwhm_(~self.outliers)', 'norm');
        end
        
        function fit = fit_fr_peak(self)
            fr_max = self.fr_at_peak;
            fit = fitdist(fr_max(~self.outliers)', 'norm');
        end
        
        function [V, p, v0, phi, speed] = regress(self, field)
            if nargin < 2 || isempty(field), field = 'onsets'; end
            y = self.(field)(~self.outliers)';
            
            X = [ones(size(y(:))) self.position(~self.outliers, :)];
            [b, ~, ~, ~, stats] = regress(y(:), X);
            p = stats(3);
            v0 = b(1);
            V = pinv(b(2:3));
            phi = atan2(V(2), V(1));
            speed = norm(V);
            
        end
        function fit = wave_fit(self)
            
            [V_, p, v0] = self.regress;
            fit.V = V_;
            fit.speed = norm(V_) * .4;  % in mm / s  
            fit.phi = atan2(V_(2), V_(1));
            fit.p = p;
            fit.v0 = v0;
        end
        function s_dur = seizure_dur(self)
            p = strsplit(self.name, '_');
            if strcmpi(p{1}, 'c7') 
                s_dur = 34;
            else
                load([p{1} filesep ...
                    self.name '_Neuroport_10_10.mat'], ...
                    'Duration');
                s_dur = Duration;
            end
        end
        function orel = onset_rel(self)
            dur = self.seizure_dur;
            orel = self.center./dur;
        end
        function acf = autocorr(self, winsize, offset)
            % moving peak temporal autocorr. 
            %     <window> in ms (default: 1000)
            %     <offset> in ms (default: 100)
            if nargin < 2, winsize = 1e3; end 
            if nargin < 3, offset = 100; end
%             fr_ = nanmean(self.fr(:, ~self.outliers), 2);
            fr_ = self.fr(:, ~self.outliers);
            winsize = winsize + offset;
            [Nt, Nch] = size(fr_);
            ac = nan(Nt, Nch);
            for ii = 1:100:Nt - winsize
                for jj = 1:Nch
                    inds = ii - 1 + (1:winsize);
                    [aa, lags] = xcorr(fr_(inds, jj), 'normalize');
                    ac(ii, jj) = max(aa(lags > offset));
                end
            end
            acf.trace = nanmean(ac, 2);
            acf.median_late = nanmedian(ac(self.time > self.range(2) & self.time < self.seizure_dur));
            acf.median_iw = nanmedian(ac(self.time > self.range(1) & self.time < self.range(2)));
            acf.median_all = nanmedian(ac(self.time > 0 & self.time < self.seizure_dur));
            
        end
        function ff = fano(self)
            ei = self.mea.event_inds;
            k = self.mea.SRO * 1;  % use a 1 second window
            sig2 = movvar(single(full(ei)), k);
            mu = movmean(single(full(ei)), k);
            fano_ = sig2 ./ mu;
            fano_(:, self.outliers) = [];
            fno = decimate( ...
                fillmissing(double(nanmean(fano_, 2)), 'nearest'), ...
                self.mea.SRO ./ self.mea.SamplingRate);
            ff.trace = fno;
            ff.min = median(min(fano_(...
                self.mea.AllTime > self.range(1) & ...
                self.mea.AllTime < self.range(2), :)));
            ff.median_late = nanmedian(fno(self.time > self.range(2) & self.time < self.seizure_dur));
            ff.median_iw = nanmedian(fno(self.time > self.range(1) & self.time < self.range(2)));
            ff.median_all = nanmedian(fno(self.time > 0 & self.time < self.seizure_dur));
        end
        function ff = fano2(self)
            fr_ = self.fr(:, ~self.outliers);
            k = self.mea.SamplingRate * 1;  % use a 1 second window
            fno = movvar(fr_, k) ./ movmean(fr_, k);
            ff.trace = nanmean(fno, 2);
            ff.max = median(max(fno(self.time > self.range(1) & self.time < self.range(2), :)));
            ff.median_late = nanmedian(fno(self.time > self.range(2) & self.time < self.seizure_dur));
            ff.median_iw = nanmedian(fno(self.time > self.range(1) & self.time < self.range(2)));
            ff.median_all = nanmedian(fno(self.time > 0 & self.time < self.seizure_dur));
        end
        function ff = vmr(self)
            frS = self.fr_smooth(:, ~self.outliers); 
            fno = nanvar(frS, [], 2) ./ nanmean(frS, 2);
            ff.trace = fno;
            ff.median_late = nanmedian(fno(self.time > self.range(2) & self.time < self.seizure_dur));
            ff.median_iw = nanmedian(fno(self.time > self.range(1) & self.time < self.range(2)));
            ff.median_all = nanmedian(fno(self.time > 0 & self.time < self.seizure_dur));
        end
        function ctr = crossing_time_rel(self)
            ctr = diff(self.range) / self.seizure_dur;
        end
        function er = ending_rel(self)
            er = self.range(2) / self.seizure_dur;
        end
        function fr = fr_median_iw(self)
            frS = nanmean(self.fr_smooth(:, ~self.outliers), 2);
            mask = self.time > self.range(1) & self.time < self.range(2);
            fr = nanmedian(frS(mask));
        end
        function fr = fr_median_all(self)
            frS = nanmean(self.fr_smooth(:, ~self.outliers), 2);
            mask = self.time > 0 & self.time < self.seizure_dur;
            fr = nanmedian(frS(mask));
        end
        function fr = fr_median_late(self)
            frS = nanmean(self.fr_smooth(:, ~self.outliers), 2);
            mask = self.time > self.range(2) & self.time < self.seizure_dur;
            fr = nanmedian(frS(mask));
        end
        function pp = fr_peak_prom(self)
            fr_ = nanmean(self.fr_smooth(:, ~self.outliers), 2);
            frS = smoothdata(fr_, 'gaussian', 5, 'samplepoints', self.time);
            [~, ~, ~,pp] = findpeaks(frS, 'NPeaks', 1, 'SortStr', 'descend');
        end
        
        function mi = moransI(self)
            ons = self.onsets(~self.outliers);
            pos = self.position(~self.outliers, :);
            mi = moran_index(ons, pos);
        end
        function mi = moransI_fr(self)
            fr_ = self.fr(:, ~self.outliers);
            N = size(fr_, 1);
            trace = nan(N, 1);
            pos = self.position(~self.outliers, :);
            for ii = 1:N, trace(ii) = moran_index(fr_(ii, :), pos); end
            mi.trace = trace;
            mi.median_late = nanmedian(trace(self.time > self.range(2) & self.time < self.seizure_dur));
            mi.median_iw = nanmedian(trace(self.time > self.range(1) & self.time < self.range(2)));
            mi.median_all = nanmedian(trace(self.time > 0 & self.time < self.seizure_dur));
        end
        function mi = moransI_frS(self)
            frS = self.fr_smooth(:, ~self.outliers);
            N = size(frS, 1);
            trace = nan(N, 1);
            pos = self.position(~self.outliers, :);
            for ii = 1:N, trace(ii) = moran_index(frS(ii, :), pos); end
            mi.trace = trace;
            mi.median_late = nanmedian(trace(self.time > self.range(2) & self.time < self.seizure_dur));
            mi.median_iw = nanmedian(trace(self.time > self.range(1) & self.time < self.range(2)));
            mi.median_all = nanmedian(trace(self.time > 0 & self.time < self.seizure_dur));
        end
        
        
        function [fit, Vx, Vy, toa] = wave_fit_alt(self)

            pos = self.position;
            OTL = self.outliers(:);
            Vx = nan(10);
            Vy = nan(10);
            toa = nan(10);
            for ii = 2:9
                for jj = 2:9
                    mask1 = pos(:, 1) >= ii - 1 & pos(:, 1) <= ii + 1;
                    mask2 = pos(:, 2) >= jj - 1 & pos(:, 2) <= jj + 1;
                    mask = ~(mask1 & mask2) | OTL;
                    self.ManualOuts = mask;
                    if sum(~mask) < 5, continue; end
                    [V_, p, v0] = self.regress;
%                     fprintf('(%d, %d): %.04g\n', ii, jj, p)
                    if ~(p < .05), continue; end
                    Vx(ii, jj) = V_(1);
                    Vy(ii, jj) = V_(2);
                    toa(ii, jj) = v0;
                end
            end
            mag = vecnorm([Vx(:) Vy(:)], 2, 2);
            fit.speed = nanmean(mag) * .4;  % in mm / s  
            fit.V = nanmean([Vx(:) Vy(:)]);
            fit.Vn = nanmean([Vx(:)./mag Vy(:)./mag]);
            fit.phi = atan2(fit.V(2), fit.V(1));
            fit.dir_coeff = vecnorm(fit.Vn);
            self.ManualOuts = [];
            
            
            
%             fit.V = V_;
%             fit.phi = atan2(V_(2), V_(1)) / pi * 180;
%             fit.p = p;
%             fit.v0 = v0;
            
        end
        function [pos, temp, zz, mag, div, crl] = find_smooth_flow(self, thresh, h)
            % [t2, zz, magS, pos] = self.find_smooth_flow(thresh=[.7 .7 .7]);
            % Input: thresh = [mag div curl] thresholds
            % Isolates spatial clusters of electrodes with similar flow
            % (2D gradient), or high levels of divergence/curl.
            
            if nargin < 2 || isempty(thresh), thresh = [.5 .5 .5]; end
            if nargin < 3, h = []; end
            mag_thresh = thresh(1);
            div_thresh = thresh(2);
            crl_thresh = thresh(3);
            
            
            % compute gradient
            [fy, fx, zz] = self.gradient;
            
            
            mag = sqrt(fx.^2 + fy.^2);  % ... and norm of gradient
            div = divergence(fy, fx);
            crl = curl(fy, fx);
            mask1 = mag < mag_thresh;
            mask2 = abs(div) < div_thresh;
            mask3 = abs(crl) < crl_thresh;
            temp = nan(10);
            temp(self.locs(~self.outliers)) = self.onsets(~self.outliers);
            temp(mask1 & mask2 & mask3) = nan;
            loners = conv2(~isnan(temp), ones(3), 'same') < 3;  % exclude electrodes with fewer than two neighbor
            temp(loners) = nan;
            
            if nargout > 0
                [px, py] = find(~isnan(temp));
                pos = [px(:) py(:)];
            end
            
            if nargout == 0 && isempty(h), h = figure; end
            if ~isempty(h)
                set(0, 'currentfigure', h);
                imagesc(temp); colorbar;
                title({strrep(self.name, '_', ' '); ...
                    ['wave ' num2str(self.wave)]});
            end
            
        end
        function [fy, fx, zz] = gradient(self)
            temp = nan(10);
            temp(self.locs(~self.outliers)) = self.onsets(~self.outliers);

            [x, y] = find(~isnan(temp));
            z = normalize(temp(~isnan(temp)), 'center');
            

            F = scatteredInterpolant( ...
                x, y, z, 'natural');
            [xx, yy] = ndgrid(0:11, 0:11);
            zz = F(xx, yy);
            
            [fy, fx] = gradient(zz);  % compute gradient
            zz = zz(2:11, 2:11);
            fy = conv2(fy, ones(3)/9, 'valid');
            fx = conv2(fx, ones(3)/9, 'valid');
            mask = isnan(temp);
            fy(mask) = nan;
            fx(mask) = nan;
            
            
        end
        function sp = speed_(self)
            fit = self.wave_fit_alt;
            sp = fit.speed;
        end
        function show = get.show(self)
            show = sum(~self.outliers) >= self.MinElectrodes;  % ...
%                 && self.center/self.mea.Duration < .5;
        end

        
        
        function ax = plot2D(self, type, ax)
            if nargin < 3, ax = gca; end
            if nargin < 2 || isempty(type), type = 'onsets'; end
            if isa(type, 'matlab.graphics.axis.Axes'), ax = type; type = 'onsets'; end
            type = validatestring(type, { ...
                'onsets', 'power', 'firingrate', 'duration' ...
                });
            
            switch type
                case 'onsets'
                    axis(ax, 'square')
                    set(ax, 'nextplot', 'replacechildren', 'units', 'points');
                    gs = self.GridSize;
                    pt_width = min(ax.Position([3 4])) / max(gs);
                    
                    
                    
                    fit = self.wave_fit;
                    ons = self.onsets(~self.outliers);
                    pos = self.position(~self.outliers, :);
                    sz = self.fr_at_peak(~self.outliers);
                    szR = rescale(sz, (pt_width / 4).^2, (1.1*pt_width).^2);
                    scatter(ax, pos(:, 1), pos(:, 2), szR, ons, 'filled', 's');
                    
                    xlim([0 gs(1)+1]);
                    ylim([0 gs(2)+1]);
                    colorbar(ax);
                    title(ax, {strrep(self.name, '_', ' '); 'IW onset time (s)'});
                    
                    
                    % Add direction arrow
                    V_ = .8 * gs(1) * fit.V./vecnorm(fit.V);
                    x = V_/2;
                    hold on
                    quiver(ax, (gs(1) + 1)/2 - x(1), (gs(2) + 1)/2 - x(2), ...
                        V_(1), V_(2), 0, ...
                        'color', [0 0 0], 'linewidth', 2)
                    hold off
                    xstring = strsplit(sprintf('Speed (mm/s): %0.3f, phi: %0.1f, p<5e%d', ...
                        fit.speed, atan2(V_(2), V_(1)) / pi * 180, ceil(log10(fit.p/5))), ',')';
                    xlabel(ax, xstring); 
                    
                case 'power'
                    % Peak power during IW
                    data = self.peak_pwr;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar(ax);
                    title(ax, {strrep(self.name, '_', ' '); 'Peak power'});
                    
                    axis(ax, 'square')
                    
                case 'firingrate'
                    % Peak power during IW
                    data = self.fr_at_peak;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Firing rate at peak (spikes/s)'});
                    
                    axis(ax, 'square')
                    
                case 'duration'
                    % Peak power during IW
                    
                    data = self.durs;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Duration of IW [s]'});
                    
                    axis(ax, 'square')

            end
        end
        function ksdensity(self, type)
            if nargin < 2, type = 'onsets'; end
            type = validatestring(type, {'onsets', 'durs', 'power', 'freq'});
            
            
            switch type
                case 'onsets'
                    ksdensity(self.all_locs, ...
                        'NumPoints', 2*ceil(range(self.time)), ...
                        'bandwidth', 2); %#ok<CPROPLC>
                    hold on
                    ksdensity(self.onsets(~self.outliers), ...
                        'NumPoints', 2*ceil(range(self.time)), ...
                        'bandwidth', 2); %#ok<CPROPLC>
                    hold off
                    title({strrep(self.name, '_', ' '); ...
                        'Onset times'});
                    xlabel('Time (s)')
                    
                case 'durs'
                    ksdensity(self.all_durs, ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection');
                    hold on
                    ksdensity(self.durs(~self.outliers), ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection');
                    hold off
                    title({strrep(self.name, '_', ' '); ...
                        'Durations'});
                    xlabel('Time (s)')
                    
                case 'power'
                    pwr = self.peak_pwr + eps;
                    ksdensity(pwr, ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection');
                    hold on
                    ksdensity(pwr(~self.outliers), ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection')
                    hold off
                    title({strrep(self.name, '_', ' '); ...
                        'Peak power'});
                    xlabel('Power')
                    
                case 'freq'
                    frq = self.fr_frq;
                    ksdensity(frq, ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection');
                    hold on
                    ksdensity(frq(~self.outliers), ...
                        'support', 'positive', ...
                        'boundarycorrection', 'reflection');
                    hold off
                    title({strrep(self.name, '_', ' '); ...
                        'Peak frequencies'});
                    xlabel('Freq (Hz)')

            end
            
            legend('All', 'Recruited');
            ylabel('PDF')

        end
        function plot(self, style, ax, show_onsets)
            % self.plot(style='raster');  
            % Valid styles are {'raster', 'lowpass'}. 
            %   lowpass: full (normalized) low pass trace on each channel.
            %   raster: shows suprathreshold (self.MinPeakHeight) intervals
            %       on each channel with red dot indicating wave onset time
            
            if nargin < 2 || isempty(style), style = 'raster'; end  
            if nargin < 3 || isempty(ax), ax = gca; end
            if nargin < 4, show_onsets = false; end
            
            style = validatestring(style, ...
                {'raster', 'lowpass', 'raster_fr', 'fr', 'fr_pow'});
            switch style
                case 'fr_pow'
                    fr_mn = nanmean(self.fr, 2);
                    plo_mn = nanmean(self.p_lo, 2);
                    yyaxis(ax, 'left')
                    plot(ax, self.time, fr_mn, ...
                        'linewidth', .5, 'color', .8 * [1 1 1])
                    ylim(ax, [0 1.1*max(fr_mn)]);
                    ylabel(ax, 'Firing rate (spikes/s)');
                    yyaxis(ax, 'right')
                    plot(ax, self.t, plo_mn, ...
                        'linewidth', 2, 'color', .15 * [1 1 1])
                    ylim(ax, [0 1.1*max(plo_mn)]);
                    ylabel(ax, {'Power [0-2] Hz'; '(normalized)'})
                    
                    for ii = 1:self.num_waves
                        self.wave = ii;
                        if sum(isfinite(self.onsets)) < 10, continue; end
                        t0 = nanmedian(self.onsets);
%                         wave_range = ...
%                             [self.onsets(:) - self.durs(:)/2, ...
%                             self.onsets(:) + self.durs(:)/2];
%                         x_range = [min(wave_range(~self.outliers, 1)), ...
%                             max(wave_range(~self.outliers, 2))];
                        x_range = self.range;
                        y_range = ylim(ax);
                        xline(ax, t0);
                        
                        [v1, v2] = ndgrid(x_range, y_range);
                        pp = patch(ax, 'vertices', [v1(:) v2(:)], ...
                            'faces', [1 2 4 3], ...
                            'facecolor', .5*[1 1 1], ...
                            'facealpha', .5, ...
                            'linestyle', 'none'); %#ok<NASGU>
                    end

                    axis(ax, 'tight')
                    title(ax, {strrep(self.name, '_', ' ')});
                case 'fr'
                    
                    fr_mn = mean(self.fr, 2);
                    fr_mn1 = mean(self.fr_smooth, 2);
                    ff = self.vmr;
                    ff = ff.trace;

                    plot(ax, self.time, rescale(fr_mn, 0, max(fr_mn1(:))), ...
                        'linewidth', .5, 'color', .8 * [1 1 1], 'tag', 'fr')
                    ylabel(ax, {'Firing rate'; '(spikes/s)'});
                    hold(ax, 'on')
                    plot(ax, self.time, fr_mn1, ...
                        'linewidth', 2, 'color', .15 * [1 1 1], 'tag', 'frS');
                    
                    plot(ax, self.time, ff, ...
                        'linewidth', 2, 'color', .4 * [1 1 1], 'tag', 'VMR');
                    hold(ax, 'off')
                    
                    ylabel(ax, {'Firing rate'; '(spikes/s)'})
                    
                    for ii = 1:self.num_waves
                        self.wave = ii;
                        if ~self.show, continue; end
                        t0 = self.center;
                        x_range = self.range;
                        y_range = ylim(ax);
                        xline(ax, t0);
                        
                        [v1, v2] = ndgrid(x_range, y_range);
                        pp = patch(ax, 'vertices', [v1(:) v2(:)], ...
                            'faces', [1 2 4 3], ...
                            'facecolor', .5*[1 1 1], ...
                            'facealpha', .5, ...
                            'linestyle', 'none'); %#ok<NASGU>
                    end

                    axis(ax, 'tight')
                    dat = [fr_mn1, ff];
                    dat = dat(self.time > 0 & self.time < self.seizure_dur, :);
                    ylim(ax, [0 1.1*max(dat(:))])
                    title(ax, {strrep(self.name, '_', ' ')});
                case 'raster_fr'
                    frS = self.fr_smooth;
                    frS(self.time < 0 | self.time > self.seizure_dur, :) = nan;
                    frN = normalize(frS);
%                     [n_t, n_ch] = find(frN > self.MinPeakHeight & frS > self.MinPeakFr);
%                     plot(ax, self.time(n_t), n_ch, '.');
%                     hold on
%                     plot(self.onsets, 1:self.nch, 'r.');
%                     hold off
                    mat = frN > self.MinPeakHeight & frS > self.MinPeakFr;
                    [ons, so] = sort(self.onsets);
                    mat = mat(:, so);
%                     mat = frS > self.MinPeakFr;
%                     mat = false(size(frS));
%                     tt = self.time;
%                     for ii = 1:size(mat, 2)
%                         bounds = self.onsets(ii) + [-.5 .5] * self.durs(ii);
%                         mask = tt >= bounds(1) & tt <= bounds(2);
%                         mat(mask, ii) = true;
%                     end
                    imagesc(ax, self.time, 1:numel(self.onsets), mat');
                    if show_onsets
                        hold on
                        plot(ons,  1:self.nch, 'r.');
                        hold off;
                    end
                    colormap(ax, 1 - gray(2)); axis(ax, 'xy')
                    axis(ax, 'tight');
                    xlim(ax, [self.time(1) self.time(end)])
                    title(ax, {strrep(self.name, '_', ' '); ...
                        sprintf('Mean duration = %0.3f s', nanmean(self.durs))});
                    ylabel(ax, 'Channel')
                case 'raster'
%                     pwr = self.p_lo;
%                     inds = interp1(self.time, 1:numel(self.time), self.t, 'nearest');
%                     frS = self.fr_smooth(inds, :);
%                     [n_t, n_ch] = find(pwr > self.MinPeakHeight & frS > self.MinPeakFr);
%                     plot(ax, self.t(n_t), n_ch, '.');
%                     mat = pwr > self.MinPeakHeight & frS > self.MinPeakFr;
                    frS = self.fr_smooth;
                    mat = frS > self.MinPeakFr & normalize(frS) > self.MinPeakHeight;
                    imagesc(ax, self.time, 1:numel(self.onsets), mat');
                    if show_onsets
                        hold on
                        plot(self.onsets,  1:self.nch, 'r.');
                        hold off;
                    end
                    colormap(ax, 1-gray(2)); axis(ax, 'xy')
%                     hold on
%                     plot(ax, self.onsets, 1:self.nch, 'r.');
%                     hold off
                    axis(ax, 'tight');
                    xlim(ax, [self.time(1) self.time(end)])
                    title(ax, {strrep(self.name, '_', ' '); ...
                        sprintf('Mean duration = %0.3f s', nanmean(self.durs))});
                    ylabel(ax, 'Channel')
                case 'lowpass'
                    switch self.Method
                        case 'pwr'
                            plot(ax, self.t, self.p_lo/4 + (1:self.nch), 'k') 
                        case 'fr'
                            plot(ax, self.time, normalize(self.fr_smooth)/4 + (1:self.nch), 'k') 
                    end
                    axis(ax, 'tight')
                    title(ax, strrep(self.name, '_', ' '));
                    ylabel(ax, 'Channel')
            end
            xlabel(ax, 'Time (s)');
        end
        function reset(self)
%             [self.onsets, self.durs, self.outliers] = deal([]);
        end
        
        function h = inspect(self)
            h = figure;
            self.plot;  % *** visualize: iw candidate times ***
        end
        
        function [stats, h] = get_stats(self, plt)
            % [stats, h] = get_stats(self, plt=false)
            [stats, h] = deal([]);
            if nargin < 2, plt = false; end
            if sum(~self.outliers) < 5, disp('Skipping: fewer than 5 data points.'); return; end
            
%             self.reset;


            % Get all the stats that used to be in iw_info... (combining
            % these files...)
            stats = struct( ...
                'wave_num', self.wave, ...
                'center', self.center, ...  
                'range', self.range, ...
                'phi', self.wave_fit.phi, ...
                'phi_alt', self.wave_fit_alt.phi, ...
                'phi_pval', self.wave_fit.p, ...
                'main_wave', self.main_wave, ...
                'template', self.iw_templates);  % this will result in a copy of templates being saved for each wave, but I don't think it should be an issue since there aren't that many waves and this should be small
            
            for ff = {'prop_sp', 'moransI', 'range', ...
                    'wave_fit', 'fit_fr_peak', 'fit_iw_fwhm', 'main_wave', ...
                    'wave_fit_alt', 'crossing_time_rel', 'ending_rel', ...
                    'fr_median_iw', 'fr_median_all', 'fr_median_late', 'fr_fwhm'}
                    
                stats.(ff{:}) = self.(ff{:});
            end
            
            
            stats.onset = self.center;  % this is a duplicate of stats.center, but I don't know what calls things this way so leave it
            
            
            fano = self.fano;
            for fn = fieldnames(fano)'
                if numel(fano.(fn{:})) > 2, continue; end
                stats.(['fano_' fn{:}]) = fano.(fn{:});
            end
            
            fano2 = self.fano2;
            for fn = fieldnames(fano2)'
                if numel(fano2.(fn{:})) > 2, continue; end
                stats.(['fano_neuro_' fn{:}]) = fano2.(fn{:});
            end
            
            acf = self.autocorr;
            for fn = fieldnames(acf)'
                if numel(acf.(fn{:})) > 2, continue; end
                stats.(['autocorr_' fn{:}]) = acf.(fn{:});
            end
            
            vmr = self.vmr;
            for fn = fieldnames(vmr)'
                if numel(vmr.(fn{:})) > 2, continue; end
                stats.(['vmr_' fn{:}]) = vmr.(fn{:});
            end
            
            MI = self.moransI_fr;
            for fn = fieldnames(MI)'
                if numel(MI.(fn{:})) > 2, continue; end
                stats.(['moransI_' fn{:}]) = MI.(fn{:});
            end
            
            MI = self.moransI_frS;
            for fn = fieldnames(MI)'
                if numel(MI.(fn{:})) > 2, continue; end
                stats.(['moransI_sm_' fn{:}]) = MI.(fn{:});
            end

                        
            
            if plt
                clear h
                close all
                h(1) = figure();
                self.plot;
                h(2) = figure;
                self.plot2D;
                fc = 3;
                for type = {'onsets', 'iw_fwhm', 'power', 'freq'}
                    h(fc) = figure;
                    self.ksdensity(type{:});
                    fc = fc + 1;
                end


                h(fc) = figure;
                if PS, self.find_smooth_flow([], h(fc)); end
                fc = fc + 1;

                for type = {'fr_pow', 'fr', 'raster_fr'}
                    h(fc) = figure();
                    self.plot(type{:});
                    fc = fc + 1;
                end
            
            
            end
            
            display(stats.prop_sp, 'prop speed')
            display(stats.wave_fit, 'wave_fit')
            
            
        end

    end
    
    
    methods (Static)
        
        function outliers = outliers_by_clust(onsets, pos, ~)
        % outliers = outliers_by_clust(onsets, pos)
        % Finds outliers using time and position. Aborts if the cluster size goes
        % below 10. Clusters must be contiguous with neighboring onset times
        % differing by no more that 1/2 sd of all times.
        %
        % Considered using fr_at_peak as another measure of distance, but figured
        % best not for now so we don't wash out this information.

            while sum(isfinite(onsets)) >= 10
                X = [ ...
                    pos, ...  
                    normalize(onsets(:))]; %, ...
        %             normalize(fr(:))];
                mask = isfinite(onsets);
                tt = clusterdata(X(mask, :), ...  
                    'distance', 'euclidean', ...
                    'criterion', 'distance', ...
                    'linkage', 'single', ...  % Nearest neighbor
                    'cutoff', sqrt(2 + .5));  % cluster must be connected and values differ by less that .5 sd
                counts = histcounts(categorical(tt));
                [N, g] = max(counts);
                if numel(counts) == 1 || N < 10
                    break
                end
                ons_temp = onsets(mask);
                ons_temp(tt ~= g) = nan;
                onsets(mask) = ons_temp;
            end
            outliers = isnan(onsets);
        end
        
        function outliers = outliers_by_MAD(onsets)
            mask = isoutlier(onsets);
            while any(mask)  % remove outliers
                onsets(mask) = nan;
                mask = isoutlier(onsets);
            end
            onsets(mask) = nan;
            outliers = isnan(onsets);  % find all outliers  
        end

        function iw = loadobj(iw)
            
            iw.locs = [];  % had this save, when should have been transient            
%             iw.reset;
        end
        
        function T = add_property_to_statfile(property, sz_num, folder, stat_file)
            % T = add_property_to_statfile(property, sz_num=[], folder='iw_mats', stat_file='iw_stats')
            % Input property must be an IW property or method with a scalar
            % output value.
            if ischar(property), property = {property}; end
            if nargin < 2, sz_num = []; end
            if nargin < 3, folder = 'iw_mats'; end
            if nargin < 4, stat_file = 'iw_stats'; end
            
            
            out = SeizureInfo;
            if isempty(sz_num), sz_num = 1:height(out); end
            stat_file = matfile(stat_file, 'writable', true);
            
            for sz = sz_num
          
                fname = sprintf('%s_Seizure%d.mat', out.patient{sz}, out.seizure(sz)); 
%                 fprintf('%d: %s\n', sz, fname)
                load([folder filesep fname], 'iw');
                

                for ii = 1:iw.num_waves
                    
                    iw.wave = ii;
%                     iw.reset;
                    fieldname = [iw.name '_' num2str(iw.wave)];
                    fprintf('%d,%d: %s\n', sz, ii, fieldname)
                    
                    if sum(~iw.outliers) < 5, fprintf('Skipping %s: fewer than 5 data points.\n', fieldname); continue; end

                    stats = stat_file.(fieldname);
                    for pp = property
                        if strcmpi(pp{:}, 'vmr')
                            ff = iw.vmr;
                            for sub = fieldnames(ff)'
                                if numel(ff.(sub{:})) > 2, continue; end
                                stats.(['vmr_' sub{:}]) = ff.(sub{:});
                            end
                        else
                            stats.(pp{:}) = iw.(pp{:});
                        end
                    end
                    
                    stat_file.(fieldname) = stats;
                    
                end

            end
            
            stats = load(stat_file.Properties.Source);
            T = IW.stats2table(stats);
            
        end
        
        function T = compute_iw_properties(sz_num, PLOT, WHITEN, REMAKE_IW)
            % compute_iw_properties(sz_num=[], PLOT=false, WHITEN=false, REMAKE_IW=false)
            % Recompute IW properties and save results to <iw_mats/> and
            % <iw_stats.mat>
            if nargin < 1 || isempty(sz_num), sz_num = []; end
            if nargin < 2 || isempty(PLOT), PLOT = false; end
            if nargin < 3 || isempty(WHITEN), WHITEN = false; end
            if nargin < 4 || isempty(REMAKE_IW), REMAKE_IW = false; end
            
            out = SeizureInfo;
            if isempty(sz_num), sz_num = 1:height(out); end
            if isa(sz_num, 'MEA'), mea = sz_num; sz_num = 0; end

            folder = 'iw_mats';
            stat_folder = 'iw_stats';
            
            if WHITEN
                folder = 'iw_mats_wh'; 
                stat_folder = 'iw_stats_wh';
            end

            for sz = sz_num
                close all
                if sz == 0
                    fprintf('%d: %s\n', sz, mea.Name)
                else
                    fname = out.fname{sz};
                    mea = MEA(fname);
                    fprintf('%d: %s\n', sz, fname)
                end


                if REMAKE_IW
                    iw = IW(mea);
                    % Remove any old files for the pat/seizure
                    delete([folder filesep mea.Name '_*.mat'])
                    delete([stat_folder filesep mea.Name '_*.mat'])
                    delete(['iw_info' filesep mea.Name '.mat'])  % as of 5/5/21, this is no longer updated; allow removal of files with the next run to ensure no old data is being used
                else    
                    load([folder filesep mea.Name], 'iw');
                    iw.mea = mea;
                end

                if WHITEN, iw.mea.whiten; end
                iw.Method = 'fr';

                for ii = 1:iw.num_waves
                    
                    iw.wave = ii;
%                     iw.reset;

                    fieldname = [iw.name '_' num2str(iw.wave)];


                    [stats, h] = iw.get_stats(PLOT);
                    save([stat_folder filesep fieldname], 'stats');

%                     if ii == iw.main_wave; BVNY.create_iw_info(iw); end
                    ss = mkdir([folder filesep fieldname]);
                    if ~ss, error('Could not make directory.'); end
                    if ~isempty(h), savefig(h, [folder filesep fieldname]); end
                    for ih = 1:numel(h)  % print figs as PNGs
                        print(h(ih), [folder filesep fieldname filesep num2str(ih)], '-dpng'); 
                    end
                    

                end
                if REMAKE_IW
                    save([folder filesep mea.Name], 'iw'); 
                end

            end
            

            stats = foldercontents2struct(stat_folder);
            IW.stats2table(stats);
            if nargout > 0, T = readtable('iw_table'); end
        end
        
        function T = stats2table(stats)

            % Load from the default stats location. 
            if nargin < 1, stats = foldercontents2struct('iw_stats'); end
            
            rows = fieldnames(stats);
            
            % Only use fields that appear in all non-empty rows
            mask = cellfun(@(x) isempty(stats.(x)), rows);
            temp = cellfun(@(x) fieldnames(stats.(x)), rows(~mask), 'uni', 0); 
            [hh, c] = histcounts(categorical(cat(1, temp{:})));
            fields = c(hh == sum(~mask));
            if any(hh < sum(~mask))
                fprintf('Warning: The following fields are not found in all rows:\n');
                exc_f = sprintf('%s, ', c{hh < sum(~mask)});
                fprintf('%s\n', exc_f(1:end-2))
            end
            
            for ii = numel(rows):-1:1
                ff = rows{ii};
                dat = stats.(ff);
                if isempty(dat) || ~isfield(dat, 'prop_sp'), continue; end
                fieldinf = strsplit(ff, {'Seizure', '_'});
                row.patient = fieldinf{1};
                row.seizure = str2double(fieldinf{2});
                row.wave_num = str2double(fieldinf{3});

%                 row.onset = dat.onset;
%                 row.prop_sp = dat.prop_sp;
%                 row.t_bounds = dat.t_bounds;
                
                % Get most fields directly; manage a few individually
                % later.
                for fd = fields
                    if contains(fd{:}, {'wave_fit', 'fit_'}), continue; end
                    row.(fd{:}) = dat.(fd{:});
                end
                
                wf = dat.wave_fit;
                row.V = wf.V;
                row.speed = wf.speed;
                row.phi = wf.phi;
                row.p = wf.p;


                wf = dat.wave_fit_alt;
                row.V_alt = wf.V;
                row.speed_alt = wf.speed;
                row.phi_alt = wf.phi;
                row.dir_coeff_V_alt = wf.dir_coeff;
                

                load([row.patient filesep ...
                    row.patient '_Seizure' num2str(row.seizure) ...
                    '_Neuroport_10_10.mat'], ...
                    'Duration');
                row.onset_rel = dat.onset / Duration;
%                 row.ending_rel = dat.ending_rel;
                row.nchannels = numfinite(dat.fit_fr_peak.InputData.data);
%                 row.crossing_time_rel = dat.crossing_time_rel;
%                 row.nchan_smooth = dat.nchannels_sm;
                for fff = ["fr_peak" "iw_fwhm"]  % , 'ht', 'frq' (these are related to pxx method of computing IW)
                    temp_fit = dat.(sprintf('fit_%s', fff));
                    if strcmpi(fff, 'ht'), name = 'pwr'; else, name = char(fff); end
                    row.([name '_mu']) = temp_fit.mu;
                    ci = temp_fit.paramci;
                    row.([name '_mu_conf']) = temp_fit.mu - ci(1, 1);
                    row.([name '_sigma']) = temp_fit.sigma;
                    row.([name '_sigma_conf']) = temp_fit.sigma - ci(1, 2);
                end


                T(ii) = row;
            end
            
            iw_table = struct2table(T, 'asarray', 1);
            iw_table(cellfun(@isempty, iw_table.patient), :) = [];
            writetable(iw_table);
            if nargout > 0, T = readtable('iw_table.txt'); end

        end
        
        

    end
    
    
end


%% Locals


