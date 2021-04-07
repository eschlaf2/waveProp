classdef IW < handle
    
    properties
        
        name
        fmax = 2
        wave = 1
        MinPeakHeight = 2  % Defines the threshold for candidate IW crossing times in standard deviations
        MinPeakFr = 30
        W
        Method = 'fr'  % 'pwr' or 'fr'
    end
    
    properties (Transient = true)  % Change to hidden if you start using this and want to save it
        pxx = struct('W', [], 'pxx', [])
        f
        t
        
    end
    
    properties (Transient = true)
        all_pks   % Used in get_stats plots
        all_locs
        all_durs

        mea
        locs
        onsets
        durs
        outliers
        fr
        fr_smooth
        peak_pwr
        fr_at_peak
    end
    
    properties (Dependent = true)
        time
        nch
        pks
        t_mask
        V
        speed
        phi
        t_bounds
        num_waves
        p_lo
        center
        range
        show
        position
    end
    
    methods
        
        function self = IW(mea, W)
            if nargin < 2, W = 20; end
            
            self.mea = mea;
            self.name = mea.Name;
            self.W = W;
            self.mea.params.fr_window = self.W;
%             self.fr = self.mea.firing_rate / 100;  % convert to firing rate per 10 ms (this is for my intuition only)
        end
        
        
        
        
       
        function set.wave(self, value)
            if value ~= self.wave
                self.reset;
                self.wave = value;
            end
        end
        function t = get.time(self); t = self.mea.Time; end
        function nch = get.nch(self); nch = size(self.mea.Data, 2); end
        function locs = get.locs(self) 
            if isempty(self.locs), self.locs = self.mea.locs; end
            locs = self.locs;
        end
        function t = get.t(self)
            if isempty(self.t), self.pxx; end
            t = self.t;
        end
        function f = get.f(self)
            if isempty(self.f), self.pxx; end
            f = self.f;
        end
        
        
        
        function mea = get.mea(self)
            if isempty(self.mea)
                pat = strsplit(self.name, '_');
                self.mea = MEA(sprintf('%s/%s_Neuroport_10_10.mat', ...
                    pat{1}, self.name));
                self.mea.params.fr_window = self.W;
            end
            mea = self.mea;
        end
        function fr = get.fr(self)
            if isempty(self.fr) || self.W ~= self.mea.params.fr_window
                self.mea.firing_rate = [];
                self.mea.params.fr_window = self.W;
                self.fr = self.mea.firing_rate;
                self.fr_smooth = [];
            end
            fr = self.fr;
        end
        function t_mask = get.t_mask(self)
            bounds = self.t_bounds;
            t_mask = self.t > bounds(1) & self.t < bounds(2);
        end
        function pks = get.pks(self)
%             if isempty(self.pks)
                [pks, loc, dur] = self.findpeaks(); 

                mask = cellfun(@isempty, pks);
                if any(mask)
                    pks(mask) = {nan};
                    loc(mask) = {nan};
                    dur(mask) = {nan};
                end

                pks = cellfun(@(x) x(1), pks);
                self.onsets = cellfun(@(x) x(1), loc); 
                self.durs = cellfun(@(x) x(1), dur);
                
                mask = self.fr_at_peak < self.MinPeakFr;
                self.onsets(mask) = nan;
                self.durs(mask) = nan;
                
                
%                 self.pks = pks;
%                 self.onsets = onsets;
%                 self.durs = durs;
%             end
%             pks = self.pks;
 
        end
        function fr = get.fr_at_peak(self)
            inds = interp1(self.time, 1:numel(self.time), self.onsets, ...
                'nearest');

            fr = nan(size(inds));
            for ii = 1:numel(fr)
                if isnan(inds(ii)), continue; end
                fr(ii) = self.fr_smooth(inds(ii), ii);
            end
            self.fr_at_peak = fr;
            
        end
        function onsets = get.onsets(self)
            if isempty(self.onsets), self.pks; end
            onsets = self.onsets;
        end
        function durs = get.durs(self)
            if isempty(self.durs), self.pks; end
            durs = self.durs;
        end
        function outliers = get.outliers(self)
            
            if isempty(self.outliers)
                if 0 && strcmpi(self.name(1:end-1), 'CUCX4_Seizure') && self.wave == 1  
                    % I used this to make IW versions of direction
                    % analysis, but result looked similar to that where
                    % there was no electrode exclusion (apart from
                    % <BadChannels>).
                    pos = self.position;
                    outliers = pos(:, 2)' > 4 | (pos(:, 2)' > 3 & pos(:, 1)' > 7);
                else
%                     outliers = outliers_by_MAD_(self.onsets);
                    outliers = IW.outliers_by_clust(self.onsets, self.position, self.fr_at_peak);
                    if sum(~outliers) < 30
                        outliers = IW.outliers_by_MAD(self.onsets);
                    end
                end
                self.outliers = outliers;

            end
            outliers = self.outliers;
            
        end
        function [pxx_fr, f] = pxx_fr(self)
            [sxx, f] = pspectrum(self.fr, self.time);
            pxx_fr = pow2db(sxx);
        end        
        function frq = fr_frq(self)
            [pxx_, f_] = self.pxx_fr;
            [~, frq] = arrayfun( ...
                @(ich) findpeaks(pxx_(f_ < 10, ich), f_(f_ < 10), 'NPeaks', 1), ...
                1:self.nch, 'uni', 0);
            nanmask = cellfun(@isempty, frq);
            frq(nanmask) = {nan};
            frq = cat(2, frq{:});
        end
        function pp = get.peak_pwr(self)
            if isempty(self.peak_pwr)
                pp = nan(self.nch, 1);
                for ich = 1:self.nch
                    pwr_temp = squeeze(sum(self.pxx(:, :, ich)));
                    pp(ich) = max(pwr_temp);
                end
                self.peak_pwr = pp;
            end
            pp = self.peak_pwr;
        end
        function V = get.V(self)
            [V, p] = self.regress;
%             if p >= 0.05, V = [nan nan]; end
        end
        function speed = get.speed(self)
            fit = self.wave_fit;
            speed = fit.speed;
        end
        function phi = get.phi(self)
            phi = atan2(self.V(2), self.V(1));
        end
        function bounds = get.t_bounds(self)
            bounds = self.wave_bounds.(self.name);
            bounds = bounds(self.wave, :);
        end
        function num_waves = get.num_waves(self)
            all_bounds = self.wave_bounds.(self.name);
            num_waves = size(all_bounds, 1);
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
            [p1, p2] = ind2sub([10, 10], self.locs);
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

        
        
        function pxx = get.pxx(self)
            % Computes the spectrogram of the firing rate on each channel.
            % Spectrogram parameters: T = 10e3 (10 s); overlap = 9900 (9.9
            % s); freq = 0:.01:self.fmax; 

            if ~isstruct(self.pxx)
                self.pxx = struct('pxx', self.pxx, 'W', self.W);
            end
           
            if isempty(self.pxx.pxx) || self.W ~= self.pxx.W
                % recompute pxx if W changes
                disp('Computing pxx ...')
                self.pxx.W = self.W;
                self.f = []; self.t = [];
                sgram = @(x) spectrogram(x, 10e3, 9900, 0:.01:self.fmax, 1e3);
                [~, self.f, tt, pxx(:, :, self.nch)] = sgram(self.fr(:, self.nch));
                self.t = tt + self.time(1);

                for ich = 1:self.nch-1
                    [~, ~, ~, pxx(:, :, ich)] = sgram(self.fr(:, ich));
                end
                self.pxx.pxx = pxx;
                disp('Done.')
            end
            pxx = self.pxx.pxx;
        end
        function frS = get.fr_smooth(self)
            % Smooth firing rate using a 1 second moving mean
            if isempty(self.fr_smooth)
%                 frS = smoothdata(self.fr, 'gaussian', 1, ...  % use a 1 second moving mean smoother
%                     'SamplePoints', self.time);
                frS = smoothdata(self.fr, 'movmean', 1, ...
                    'SamplePoints', self.time);
                self.fr_smooth = frS;
            end
            frS = self.fr_smooth;
        end
        function p_lo = get.p_lo(self)
            % normalize power (pxx) by standard deviation to highlight temporal
            % increases and then sum over all frequencies to generate a low pass signal
            
            if 1  % method 1 (using power)
                pxxN = self.pxx ./ nanstd(self.pxx, [], 2);
                p_lo = squeeze(nanmean(pxxN)); 
            else  % method 2 (gaussian smoothing window)
                t_inds = interp1(self.t, 1:numel(self.t), self.time, ...
                    'nearest', 'extrap'); %#ok<UNRCH>
                p_lo = normalize(self.fr_smooth(t_inds, :), 'scale');
            end
                
        end
        
        
        
        function sp = stop_points(self)
            [~, all_locs, ~] = self.findpeaks();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Visual inspection for multiple waves
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            all_locs = cat(2, all_locs{:});
            [dens, xi] = ksdensity(all_locs, ...
                'numpoints', ceil(range(self.time)), 'bandwidth', 2);
            [~, sp] = findpeaks(-dens, xi);
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
            sp = range(self.onsets(~self.outliers));
        end
        function fit = fit_dur(self)
            fit = fitdist(self.durs(~self.outliers)', 'norm');
        end
        function fit = fit_ht(self)
            pwr = self.peak_pwr;
            fit = fitdist(pwr(~self.outliers), 'norm');
        end
        function fit = fit_frq(self)
            frq = self.fr_frq;
            fit = fitdist(frq(~self.outliers)', 'norm');
        end
        function fit = fit_fr_peak(self)
            fr_max = self.fr_at_peak;
            fit = fitdist(fr_max(~self.outliers)', 'norm');
        end
        function fit = fit_fr_mean(self)
            fr_mn = self.mean_fr;
            fit = fitdist(fr_mn(~self.outliers)', 'norm');
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
            % for some stupid reason I put the angle in degrees here...
            [V_, p, v0] = self.regress;
            fit.V = V_;
            fit.speed = norm(V_) * .4;  % in mm / s  
            fit.phi = atan2(V_(2), V_(1)) / pi * 180;
            fit.p = p;
            fit.v0 = v0;
        end
        function s_dur = seizure_dur(self)
            p = strsplit(self.name, '_');
            if strcmpi(p{1}, 'c7'), 
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
        function fwhm = fr_fwhm(self)
            fr_ = nanmean(self.fr_smooth, 2);
            frS = smoothdata(fr_, 'gaussian', 5, 'SamplePoints', self.time);
            [~, ~, fwhm] = findpeaks(frS, self.time, 'NPeaks', 1, 'SortStr', 'descend');
            fwhm = fwhm / self.seizure_dur;
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
                    self.outliers = mask;
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
            self.outliers = OTL';
            
            
            
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
            show = sum(~self.outliers) >= 30;  % ...
%                 && self.center/self.mea.Duration < .5;
        end

        
        
        function ax = plot2D(self, type, ax)
            if nargin < 3, ax = gca; end
            if nargin < 2, type = 'onsets'; end
            type = validatestring(type, { ...
                'onsets', 'power', 'firingrate', 'duration' ...
                });
            
            switch type
                case 'onsets'
                    axis(ax, 'square')
                    set(ax, 'nextplot', 'replacechildren', 'units', 'points');
                    pt_width = min(ax.Position([3 4])) / 11;
                    
                    cmap = [1 1 1; parula];
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        self.onsets(~self.outliers);
                    
                    fit = self.wave_fit;
%                     imagesc(ax, temp);
                    ons = self.onsets(~self.outliers);
                    pos = self.position(~self.outliers, :);
                    sz = self.fr_at_peak(~self.outliers);
                    szR = rescale(sz, (pt_width / 4).^2, (1.1*pt_width).^2);
                    scatter(ax, pos(:, 1), pos(:, 2), szR, ons, 'filled');
                    xlim([0 11]);
                    ylim([0 11]);
                    colorbar(ax);
                    title(ax, {strrep(self.name, '_', ' '); 'IW onset time (s)'});
                    colormap(ax, cmap);
                    
                    
                    % Add direction arrow
                    V_ = 8 * fit.V./vecnorm(fit.V);
                    x = V_/2;
                    hold on
                    quiver(ax, 5.5 - x(1), 5.5 - x(2), V_(1), V_(2), 0, ...
                        'color', [0 0 0], 'linewidth', 2)
                    hold off
                    xlabel(ax, sprintf('Speed (mm/s): %0.3f, phi: %0.1f, p: %0.4g', ...
                        fit.speed, atan2(V_(2), V_(1)) / pi * 180, fit.p)); 
                    
                case 'power'
                    % Peak power during IW
                    cmap = [1 1 1; parula];
                    data = self.peak_pwr;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar(ax);
                    title(ax, {strrep(self.name, '_', ' '); 'Peak power'});
                    colormap(ax, cmap);
                    axis(ax, 'square')
                    
                case 'firingrate'
                    % Peak power during IW
                    cmap = [1 1 1; parula];
                    data = self.fr_at_peak;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Firing rate at peak (spikes/s)'});
                    colormap(ax, cmap);
                    axis(ax, 'square')
                    
                case 'duration'
                    % Peak power during IW
                    cmap = [1 1 1; parula];
                    data = self.durs;
                    temp = nan(10);
                    temp(self.locs(~self.outliers)) = ...
                        data(~self.outliers);
                    imagesc(ax, temp, ...
                        quantile(data(~self.outliers), [0 .97]));
                    colorbar;
                    title(ax, {strrep(self.name, '_', ' '); 'Duration of IW [s]'});
                    colormap(ax, cmap);
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
                        'bandwidth', 2);
                    hold on
                    ksdensity(self.onsets(~self.outliers), ...
                        'NumPoints', 2*ceil(range(self.time)), ...
                        'bandwidth', 2);
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
                            'linestyle', 'none');
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
                            'linestyle', 'none');
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
            [self.onsets, self.durs, self.outliers] = deal([]);
        end
        
        
        
        function h = inspect(self)
            h = figure;
            self.plot;  % *** visualize: iw candidate times ***
            display(self.stop_points, 'stop_points');
        end
        function [stats, h] = get_stats(self, plt)
            % [stats, h] = get_stats(self, plt=false)
            [stats, h] = deal([]);
            if nargin < 2, plt = false; end
            if sum(~self.outliers) < 5, disp('Skipping: fewer than 5 data points.'); return; end
            
            self.reset;
            for ff = {'prop_sp', 'fit_dur', 'moransI', ...
                    't_bounds', 'wave_fit', 'fit_fr_peak', 'fit_fr_mean', ...
                    'wave_fit_alt', 'crossing_time_rel', 'ending_rel', ...
                    'fr_median_iw', 'fr_median_all', 'fr_median_late', 'fr_fwhm'}
                    % 'fit_frq', 'fit_ht'  % these use pxx which takes a
                    % long time. Add back in later if you change to pxx
                    % method
                stats.(ff{:}) = self.(ff{:});
            end
            
            
            stats.onset = self.center;
            
            
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

            
%             try
%                 pos = self.find_smooth_flow;
%                 stats.nchannels_sm = size(pos, 1);
%                 PS = true;
%             catch ME
%                 disp('find_smooth_flow failed')
%                 disp(ME)
%                 stats.nchannels_sm = 0;
%                 PS = false;
%             end
            
            
            if plt
                clear h
                FW = 2.73;
                close all
                h(1) = figure();
                self.plot;
                h(2) = figure;
                self.plot2D;
                fc = 3;
                for type = {'onsets', 'durs', 'power', 'freq'}
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
            display(stats.fit_dur, 'durations')
%             display(stats.fit_ht, 'power')
%             display(stats.fit_frq, 'frequency')
            display(stats.t_bounds, 't_bounds')
            display(stats.wave_fit, 'wave_fit')
%             display(stats.nchannels_sm, 'nchan smooth')
            
            
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
            iw.wave = IW.main_wave(iw.name);
            
            iw.reset;
        end
        
        function add_property_to_statfile(property, sz_num, folder, stat_file)
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
                    iw.reset;
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

            folder = 'iw_mats';
            stat_file = matfile('iw_stats', 'writable', true);
            
            if WHITEN
                folder = 'iw_mats_wh'; 
                stat_file = matfile('iw_stats_wh', 'writable', true);
            end

            for sz = sz_num
                close all
                fname = out.fname{sz};
                mea = MEA(fname);
                fprintf('%d: %s\n', sz, fname)


                if REMAKE_IW
                    iw = IW(mea);
                else    
                    load([folder filesep mea.Name], 'iw');
                    iw.mea = mea;
                end

                if WHITEN, iw.mea.whiten; end
                iw.Method = 'fr';

                for ii = 1:iw.num_waves
                    
                    iw.wave = ii;
                    iw.reset;

                    fieldname = [iw.name '_' num2str(iw.wave)];

%                     if sum(~iw.outliers) < 5
%                         fprintf('Not enough data. Skipping %d: %s\n', sz, fieldname);
%                         continue
%                     end
                    [stats, h] = iw.get_stats(PLOT);
                    stat_file.(fieldname) = stats;
                    if ii == iw.main_wave(iw.name); BVNY.create_iw_info(iw); end
                    ss = mkdir([folder filesep fieldname]);
                    if ~isempty(h), savefig(h, [folder filesep fieldname]); end
                    for ih = 1:numel(h)  % print figs as PNGs
                        print(h(ih), [folder filesep fieldname filesep num2str(ih)], '-dpng'); 
                    end
                    

                end
                if REMAKE_IW, 
                    save([folder filesep mea.Name], 'iw'); 
                end

            end
            
            stats = load(stat_file.Properties.Source);
            IW.stats2table(stats);
            if nargout > 0, T = readtable('iw_table'); end
        end
        
        function T = stats2table(stats)

            rows = fieldnames(stats);
            
            % Only use fields that appear in all non-empty rows
            mask = cellfun(@(x) isempty(stats.(x)), rows);
            temp = cellfun(@(x) fieldnames(stats.(x)), rows(~mask), 'uni', 0); 
            [hh, c] = histcounts(categorical(cat(1, temp{:})));
            fields = c(hh == sum(~mask));
            if any(hh < sum(~mask))
                warning('The following fields are not found in all rows:');
                exc_f = sprintf('%s, ', c{hh < sum(~mask)});
                fprintf('%s\b\b\n', exc_f)
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
                
                fd = fieldnames(dat); 
%                 fd = fd(contains(fd, 'vmr'));
                for f_ = fd'
                    row.(f_{:}) = dat.(f_{:});
                end
                

                load([row.patient filesep ...
                    row.patient '_Seizure' num2str(row.seizure) ...
                    '_Neuroport_10_10.mat'], ...
                    'Duration');
                row.onset_rel = dat.onset / Duration;
%                 row.ending_rel = dat.ending_rel;
                row.nchannels = numfinite(dat.fit_dur.InputData.data);
%                 row.crossing_time_rel = dat.crossing_time_rel;
%                 row.nchan_smooth = dat.nchannels_sm;
                for fff = {'dur', 'fr_peak', 'fr_mean'}  % , 'ht', 'frq' (these are related to pxx method of computing IW)
                    temp_fit = dat.(['fit_' fff{:}]);
                    if strcmpi(fff{:}, 'ht'), name = 'pwr'; else, name = fff{:}; end
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
        
        function bounds = wave_bounds()
            bounds = struct( ...
                'c4_Seizure3', [1.72 15.19], ...  % [1.72 15.19]; {1.7258   15.1871   29.6839   39.0032}  c4 is all over the place
                'c4_Seizure2', [-.5 17.85], ...  % [-.5 17.85]; {-0.5039   17.8490   31.1039}
                'c4_Seizure1', [18.9 28.7], ...  % []; {8.9830   18.8553   28.7277}
                'c3_Seizure3', [-Inf 2.85; 2.85 9.15; 9.15 Inf], ...  % [-Inf 2.85; 2.85 9.15; 9.15 Inf]; {2.8500    9.1500}
                'c3_Seizure2', [-Inf 1.33; 5 12], ...  % [-Inf 1.33; 5 12]; {1.33}
                'c3_Seizure1', [0 7], ...  % []; {7.1313   16.2250}
                'CUCX2_Seizure3', [-Inf 7.57; 7.57 27.1] + 0, ...  % [-Inf 7.57; 7.57 27.1]; {7.5794   27.0841}
                'CUCX2_Seizure2', [30 Inf] - 4, ...  % []; {10.0190   21.5025   30.1152}
                'CUCX5_Seizure6', [15 32.5; 70 Inf] + 11, ...  % [15 32.5; 32.5 66.6]; {14.9040   32.4109   66.5911   69.0921}
                'CUCX5_Seizure3', [4 22; 50 Inf] + 14, ...  % [4 22; 45 54; 54 Inf]; {4.0975   22.4200   32.8900   45.1050   54.7025}
                'BW09_Seizure3', [0 14.5; 50 58], ...  % [0 14.5; 50 58];  {14.5107   30.6733   46.8360}
                'BW09_Seizure2', [0 46.7; 46.7 72; 72 94], ...  % [0 46.7; 46.7 72; 72 94]  {46.7059   72.0269   94.8975}
                'BW09_Seizure1', [13.4 28.6; 44.7 59.1; 59.1 Inf], ...  % [13.4 28.6; 44.7 59.1; 59.1 Inf]; {13.3494   28.5921   32.1787   44.7315   59.0775   Inf}
                'CUCX4_Seizure2', [6.5 51.8; 75.4 Inf] + 6, ...  % [6.5 51.8; 75.4 Inf]; {6.4592   51.8029   61.2495   69.7515   75.4194}
                'CUCX4_Seizure1', [-Inf 13.80; 13.80 64.0; 64.0 Inf] + 2, ...  % [-Inf 13.80; 13.80 64.0, 64.0 Inf]   {13.7966   31.2897   44.1793   55.2276  (possibly also 64.0?)}
                'CUCX3_Seizure6', [13.1608, 39.8139] + 3, ...  % [13.1608, 28.1532; 28.1532 39.8139]  {13.1608   28.1532   39.8139   52.3076}
                'CUCX3_Seizure5', [0, 38.59] + 5, ...  % [0, 19.81; 19.81, 38.59]; {19.8053   38.5893   46.4160   54.2427}
                'CUCX3_Seizure4', [9.36 38.3] + 5, ...  % [9.36 26.22; 26.22, 38.3]; {9.3615   26.1154   38.3000}
                'CUCX3_Seizure3', [-Inf Inf] + 3, ...  % []; {[]}
                'CUCX3_Seizure2', [-Inf Inf] + 4, ...  % []; {[]}
                'CUCX3_Seizure1', [10 Inf] + 5, ...  % []; {[]}
                'MG63_Seizure5', [0 35], ...  % [0 35; 40 60; 70 Inf]; {34.8727   61.9717}
                'MG63_Seizure4', [0 65], ...  % [0 65; 110 Inf]; {64.3583   81.4944  101.3944}
                'MG63_Seizure3', [30 60], ...  % [30 60; 80 90; 100 Inf]; {**}
                'MG49_Seizure45', [0 15; 20 40], ...  % [0 15; 20 40]; {**}
                'MG49_Seizure43', [0 15; 19 Inf], ...  % [0 10; 19 Inf]; {**}
                'MG49_Seizure36', [0 21; 20 40], ...  % [0 20; 20 40]; {**}  Note: the early wave is not as clear as in seizures 43 and 45
                'c5_Seizure3', [14.7 31.4], ...  % []; {14.7275   31.4000   46.3175   65.6225   70.8875   77.0300   92.8250}
                'c5_Seizure2', [5 27.54; 37.1 50.4], ...  % [8 30; 40 50]; {4.6417   27.5417   37.0833   50.4417   58.0750   65.7083}
                'c5_Seizure1', [4.95, 24.42; 40.9709 53.6253], ...  % []; {4.9544   24.4228   28.3165   40.9709   53.6253}
                'c7_Seizure1', [0 21.91; 21.91 34.2] ...  % [0 21.91; 21.91 34.2]; {**}
                );
        end        
        
        function out = main_wave(name)
            % Which IW wave to use to exclude channels

            dat = struct( ...
                'c4_Seizure3', 0, ...
                'c4_Seizure2', 0, ...
                'c4_Seizure1', 0, ...
                'c3_Seizure3', 0, ...
                'c3_Seizure2', 0, ...
                'c3_Seizure1', 0, ...
                'CUCX2_Seizure3', 0, ...
                'CUCX2_Seizure2', 1, ...
                'CUCX5_Seizure6', 1, ...
                'CUCX5_Seizure3', 1, ...
                'BW09_Seizure3', 1, ...  % there is a late wave that that recruits more electrodes, but it's right before the seizure ends
                'BW09_Seizure2', 1, ...
                'BW09_Seizure1', 1, ...
                'CUCX4_Seizure2', 1, ...
                'CUCX4_Seizure1', 2, ...
                'CUCX3_Seizure6', 1, ...
                'CUCX3_Seizure5', 1, ...
                'CUCX3_Seizure4', 1, ...
                'CUCX3_Seizure3', 1, ...
                'CUCX3_Seizure2', 1, ...
                'CUCX3_Seizure1', 1, ...
                'MG63_Seizure5', 1, ...
                'MG63_Seizure4', 1, ...
                'MG63_Seizure3', 1, ...
                'MG49_Seizure45', 2, ...
                'MG49_Seizure43', 2, ...
                'MG49_Seizure36', 2, ...
                'c5_Seizure3', 1, ...
                'c5_Seizure2', 1, ...
                'c5_Seizure1', 1, ...
                'c7_Seizure1', 1 ...
                );
            out = dat.(name);
        end

    end
    
    
end









%% Old

