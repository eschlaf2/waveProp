classdef fig6 < BVNY & handle

    
properties
%     Metrics = {'M', 'E', 'D1xwh'}
    Metrics = ["M", "D10"]
    Seizures
    H
    PatientByRow = false  % if true, display each row as a patient; else each column is a patient
    Smoothing = 1  % seconds
    CompiledData
    FitDiffs
    Disagree
end

properties
    Margin = [0.25 0.25];
    Padding = 0.1;
end

methods
    function F = fig6
        F.Seizures = F.SeizureInfo;
%         mask = contains(F.SeizureInfo.class, 'TLE');
%         F.Seizures = F.SeizureInfo(mask, :);
    end
    
    
    function T = show_all_sig_corr(F, sz_name, metric, pos_or_neg)
        % Show each significant wavefit/correlation pair for a given
        % seizure number and metric
        if nargin < 2 || isempty(sz_name), sz_name = 'c7_Seizure1'; end
        if nargin < 3 || isempty(metric), metric = 'M'; end
        if nargin < 4 || isempty(pos_or_neg), pos_or_neg = "pos"; end
        
        
        % Some parameters & convenience vars
        sz_name = validatestring(sz_name, string(F.SeizureInfo.name));
        metric = string(metric);
        pos_or_neg = validatestring(pos_or_neg, ["pos" "neg" "zero"]);
        sz_num = find(contains(F.SeizureInfo.name, sz_name));
        thresh = 5e-2;
        
        % Parameters for figures
        MAKE_ONE_PAGE = true;  % only show the first page of discharges
        rows = 6; cols = 8;
        panes_per_fig = rows * cols;
        
        % Get the iw info and WaveProp object
        M = WaveProp.load(F.get_file(sz_num)); M = M.(metric);
        [iw, mw] = BVNY.get_iw_info(sz_num);
        tpl = iw{mw}.template.template(mw, :, :);
        
        
        % get the correlations and directions (only show sig fits for both)
        [rho, rho_p, rho_t] = M.correlation(tpl);
        inds = M.time2inds(rho_t);
        if MAKE_ONE_PAGE
            subset = @(so) so(1:min(numel(so), panes_per_fig));
        else
            subset = @(so) so; %#ok<UNRCH>
        end
        switch pos_or_neg
            case "pos"
                ss = find(rho_p < thresh & rho > 0 & M.p(inds) < thresh);
                [~, so] = sort(rho(ss), 'descend');
                ss = ss(subset(so));
            case "neg"
                ss = find(rho_p < thresh & rho < 0 & M.p(inds) < thresh);
                [~, so] = sort(rho(ss), 'ascend');
                ss = ss(subset(so));
            case "zero"
                ss = find(rho_p > thresh & M.p(inds) < thresh);
                [~, so] = sort(abs(rho(ss)), 'ascend');
                ss = ss(so(1:min(numel(so), panes_per_fig)));  % only show one page of these
        end
        
        % Make the figures
        
        num_figs = ceil(numel(ss) / panes_per_fig);  % max (ros*cols) panes per figure
        for ii = 1:num_figs
            h(ii) = figure('name', sprintf('rho_v_dir_%s_%d', pos_or_neg, ii)); fullwidth(true); %#ok<AGROW>
            T(ii) = tiledlayout(h(ii), rows, cols); %#ok<AGROW> % If you use 'flow' here, the scatters will have be the wrong size
            title(T(ii), sprintf("%s %s", F.SeizureInfo.display_names{sz_num}, pos_or_neg));
%             uu = h.Units; h.Units = 'inches'; width = h.Position(3)/h.Position(4); h.Units = uu;
        end
        
        
        for ii = 1:numel(ss)
            which_fig = ceil(ii / panes_per_fig);
            pane = mod(ii-1, panes_per_fig) + 1;
            ax = nexttile(T(which_fig), pane); 
            M.plot2D(ax, rho_t(ss(ii))); 
            title(ax, sprintf('t = %0.2f s', rho_t(ss(ii)) - iw{mw}.center));
            xlabel(ax, sprintf('\\rho=%0.2f, p<5e%d', ...
                rho(ss(ii)), ceil(log10(rho_p(ss(ii))/5)))); 
        end
        
        % Nice legends but they cover the plots
        lgd = findobj(T, 'type', 'legend'); 
        set(lgd, 'visible', 'off'); 
        
        for ii = 1:num_figs
            tag = sprintf('%s_%s_%s_%d', F.get_fname(sz_num), metric, pos_or_neg, ii);
            F.print(h(ii), F.prefix_better(tag));
        end
    end
    
    function show_iw_all_methods(F, sz_num)
        sz = F.SeizureInfo;
        if nargin < 2 || isempty(sz_num), sz_num = 1:numel(sz.patient); end
        if isa(sz_num, 'MEA'), mea = sz_num; sz_num = 0; end
        methods = {'pwr', 'fr', 'frg'};
        r = numel(methods); c = 5;
        
        if isfinite(str2double(getenv('SGE_TASK_ID')))  % no print during testing
            print_ = @(h) {F.print(h, F.prefix_better(h.Name)); close(h)};
        else
            print_ = @(h) [];
        end
        
        S = warning; warning off
        for ss = sz_num(:)'
            if ss ~= 0, mea = MEA(ss); end
            h = figure('name', mea.Name); fullwidth(true);
            T = tiledlayout(h, r, c);
            title(T, strrep(mea.Name, '_', ' '));
            iw = mea.IW;
            iw.MaxTemplates = c;
            
            for mm = 1:length(methods)
                iw.Method = methods{mm};
                for jj = 1:min(iw.num_waves, c)
                    plotnum = sub2ind([c r], jj, mm);
                    ax = nexttile(T, plotnum);
                    iw.wave = jj;
                    iw.plot2D(ax);
                    fr = iw.fr_at_peak;
                    ylabel(sprintf('fr = %0.0f, [%d %d]', nanmean(fr), round(quantile(fr, [0 1]))));
                    title(sprintf('%s, t=%0.2f', methods{mm}, iw.center))
                    xticklabels([]);
                    yticklabels([]);
                end
            end
            print_(h);
        end
        warning(S)

    end
    function plot_clusters(F, sz_num)
        sz = F.SeizureInfo;
        if nargin < 2, sz_num = 1:numel(sz.patient); end
        if isa(sz_num, 'MEA'), mea = sz_num; sz_num = 1; end
        for ii = 1:numel(sz_num)
            pat = sz.patient{sz_num(ii)};
            seizure = sz.seizure(sz_num(ii));
            if ~exist('mea', 'var')
                mea = MEA(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', ...
                    pat, pat, seizure));
            end
            fprintf('%d: Starting %s\n', sz_num(ii), mea.Name)
            out = mea.get_IW_templates([], [], 8);
            assignin('base', 'out', out);
            assignin('base', 'mea', mea);
            
            M = WaveProp.load(mea, {'M10'});
            M = M.M10;

            h = figure('name', [M.Name '_clust']); fullwidth(true);
            M.cluster(h);
            F.print(h, F.prefix_better(h.Name));
            close(h);
            
        end
    end
    function distance_scatters(F, sz_num, method)
        % Shows the distance (correlation distance) between discharge
        % patterns and IW patterns
        PTHRESH = 5e-3;
        sz = F.SeizureInfo;
        if nargin < 3, method = ''; end  % Use default method defined in IW if not given
        if nargin < 2 || isempty(sz_num), sz_num = 1:numel(sz.patient); end
        if isa(sz_num, 'MEA'), mea = sz_num; sz_num = 0; end
        
        
        if isfinite(str2double(getenv('SGE_TASK_ID')))  % no print during testing
            print_ = @(h) {F.print(h, F.prefix_better(h.Name)); close(h)};
        else
            print_ = @(h) [];
        end
        
        
        for ii = 1:numel(sz_num)
            
            if sz_num(ii) ~= 0
                mea = MEA(sz_num(ii));
            end
            fprintf('%d: Starting %s\n', sz_num(ii), mea.Name)
            iw = mea.IW;
            if ~isempty(method), iw.Method = method; end
            out = iw.iw_templates;
            
            assignin('base', 'out', out);
            assignin('base', 'mea', mea);
            
            
            M = WaveProp.load(mea, {'M10'});
            M = M.M10;

            
            
            
            try  % tired of re-running; catch error and rerun later
            h2 = figure('name', sprintf('IW_%s_%s', out.method, mea.Name)); fullwidth(true);
            T2 = tiledlayout(h2, 5, 4);
            
            % if there are no IW, print the blank figure and skip the rest
            if isempty(out.template), nexttile(T2); title('No IW'); print_(h2); continue; end  
            
            [rho_full, ~, rho_t, rho_sig, pval_sig] = M.correlation(out.template);
            
            nexttile(T2)
            M.distance_scatter(out);
            

            mask = isfinite(rho_sig) & pval_sig(:, 1) < PTHRESH;
            title(strrep(M.Name, '_', ' '));
            
            nexttile(T2)
            histogram(pval_sig(mask, 2));
            title(sprintf('template counts, p<5e%d', log10(PTHRESH/5)))
            
            nexttile(T2)
            dd = reshape(out.template, size(out.template, 1), []);
            dd(:, all(isnan(dd))) = [];
            dd = normalize(dd, 2, 'center');
            dd = fillmissing(dd, 'constant', 0, 2);
            imagesc(squareform(pdist(dd, 'corr')), [0 2]); colorbar
            colormap(gca, make_diverging_colormap([1 0 0; 0 0 0]));
            axis square
            title('pdist(templates)')
            
            
            nexttile(T2, [3, 1])
            
            plot(rho_t - out.time', ...
                -(rho_full)/2 + (1:size(rho_full, 2)), ...
                '.', 'color', .5*[1 1 1]);
            grid on
            hold on
            if sum(mask) > 0
                gscatter(rho_t(mask) - out.time(pval_sig(mask, 2)), ...
                    -(1 - rho_sig(mask))/2 + pval_sig(mask, 2), ...   % don't forget the negative here since you put this in axis ij
                    pval_sig(mask, 2), [], [], 6, 'off')


                if 0  % Old way of plotting moving mean
                [G, gp] = findgroups(pval_sig(mask, 2));
                hc = histcounts(categorical(G));
                if any(hc == 1)
                    mask = mask & ~ismember(pval_sig(:, 2), find(hc == 1));
                    [G, gp] = findgroups(pval_sig(mask, 2));
                end
                mvmn = splitapply(@(x, t) ...
                    interp1(t, movmean(x, 10, 'SamplePoints', t), M.time', 'linear'), ...
                    1-rho_sig(mask), rho_t(mask), G);
                mvmn = fillmissing(mvmn, 'constant', 0);
                plot(M.time - out.time(gp)', -mvmn'/2 + gp', 'r', 'linewidth', 1)
                end
            end
            mvmn = movmean( ...
                fillmissing(rho_full, 'constant', 0), ...  % fill nan values with zeros
                10, 'SamplePoints', rho_t);
            plot(rho_t - out.time', -mvmn + (1:size(mvmn, 2)), 'k', 'linewidth', 2);
            hold off
            xline(0)
            ylim([0 9])
            axis ij
            title(sprintf('\\rho v. time (p<5e%d)', log10(PTHRESH/5)))
            
                        
            nexttile(T2)
            violinplot(1-rho_sig(mask), pval_sig(mask, 2));
            bp = findobj(gca, 'type', 'scatter');
            set(bp, 'MarkerFaceAlpha', 1, 'SizeData', 6);
            ylim([-1 1])
            title('\rho by template')
            
            nexttile(T2)
            violinplot(log10(pval_sig(:, 1)/5), pval_sig(:, 2));
            bp = findobj(gca, 'type', 'scatter');
            set(bp, 'MarkerFaceAlpha', 1, 'SizeData', 6);
            title('p=5eYY by template')
            
            nexttile(T2)
            vp = violinplot(rho_t(mask), pval_sig(mask, 2));
            bp = findobj(gca, 'type', 'scatter');
            set(bp, 'MarkerFaceAlpha', 1, 'SizeData', 6);
            for jj = 1:min(numel(out.time), numel(vp))
                if isempty(vp(jj).ScatterPlot), continue; end
                tt = out.time(jj);
                col = vp(jj).ViolinColor;
                yline(tt, 'color', col, 'linewidth', 2)
            end
            ylim(quantile(M.time, [0 1]))
            title('time by template');
            
            nexttile(T2)
            plot(0:100, quantile(log10(pval_sig(:, 1)/5), linspace(0, 1, 101)));
            yline(log10(PTHRESH/5));
            grid on
            title('Percentiles log10(p/5)')
            
            
            nexttile(T2)  % Leave a blank so the templates are all on the same row
            plot(0:100, quantile(abs(1-rho_sig(pval_sig(:, 1) < PTHRESH)), linspace(0, 1, 101)))
            grid on
            ylim([0 1])
            title(sprintf('Percentiles abs(\\rho), p < 5e%d', log10(PTHRESH/5)))
            
            
            nexttile(T2)
            ksdensity(1-rho_sig, linspace(-1, 1, 1000), 'bandwidth', .05);
            if sum(mask) > 1
            hold on; 
            ksdensity(1-rho_sig(mask), linspace(-1, 1, 1000), 'bandwidth', .05); 
            hold off;
            end
            legend({'all', sprintf('p<5e%d', round(log10(PTHRESH/5)))}, ...
                'location', 'best')
            title('PDF \rho')

            
            for jj = 1:size(out.template, 1)
                nexttile(T2); 
                tpl = squeeze(out.template(jj, :, :));
                fr = squeeze(out.firing_rate(jj, :, :));
                [xx, yy] = find(isfinite(tpl) );
                inds = find(isfinite(tpl) );
                scatter(xx, yy, min(fr(inds), 200) / 2, tpl(inds), 'filled')
%                 surf(fillmissing(tpl, 'nearest'), fr, 'linestyle', 'none')
%                 imagesc(squeeze(out.template(jj, :, :))); 
                if jj == iw.main_wave, str = '*'; else, str = ''; end
                title({['Template ' num2str(jj) str]; ...
                    sprintf('t=%0.2fs', nanmedian(tpl, 'all'))})
                ylabel(sprintf('fr in [%d %d]', round(quantile(fr(inds), [0 1]))));
                xlim([0 11]);
                ylim([0 11]);
                xticklabels([]);
                yticklabels([]);
                colorbar;
                axis square
            end
            
            print_(h2);
            
            
            catch ME
                warning('error in %d: %s', sz_num(ii), mea.Name)
                disp(ME)
                for line = ME.stack'
                    fprintf('In %s (line %d)\n', line.name, line.line)
                end
%                 rethrow(ME)
                print_(h2);                
            end
            
        end
%         F.print(h, F.prefix(h.Name));
    end
    
    function direction_summary_plots(F, sz_num)
        % Direction dots overlayed on hsv-colored rasters; also shows
        % dtheta/dt, directionality, discharge rate
        sz = SeizureInfo;
        if nargin < 2, sz_num = 1:numel(sz.patient); end
        
        for ii = sz_num
            pat = sz.patient{ii};
            seiz = sz.seizure(ii);
            fname = sprintf('%s_Seizure%d', pat, seiz);
            fprintf('%d: Starting %s\n', ii, fname);
            fit = WaveProp.load(fname, F.Metrics);
            [iw_info, mw] = BVNY.get_iw_info(fname);
            iw_info = iw_info{mw};
            iw_phi = iw_info.wave_fit.phi/pi*180;
            iw_t = iw_info.onset;
            
            
            for ff = F.Metrics
                W = fit.(ff{:});
                if strcmp(ff{:}, 'M10'), W.MinFinite = 10; end
                h = W.summary_plot;
                ax = findobj(h, 'Tag', 'direction_raster');

                if ~isempty(iw_info)
                    for aa = ax'
                        yline(aa, iw_phi, 'k-', 'linewidth', 2);
                    end
                
                    for aa = findobj(h, 'type', 'axes')'
                        xline(aa, iw_t, 'k-', 'linewidth', 2);
                    end
                
                pax = findobj(h, 'type', 'polar');
                hold(pax, 'on')
                polarplot(pax, [1; 1] * iw_phi, [0; 1], 'k')
                hold(pax, 'off');
                end
                
                F.print(h, F.prefix_better(sprintf('%s_%s', fname, ff{:})));
%                 close(h)
            end
        end
    end
    
    function out = transient_iw(F, interval, min_electrodes)
        % This will get all direction estimates within interval number of seconds before
        % and after the IW for each seizure and in each case return the mean direction
        % pre & post.
        
        if nargin < 2 || isempty(interval)
            fits = F.CompiledData;
        else
            fits = F.compile_fits(interval);
        end
        if nargin < 3 || isempty(min_electrodes), min_electrodes = 0; end
        
        fits = fits(fits.nchannels >= min_electrodes, :);
        
        fits.dir = angle(exp(1j*(fits.dir - fits.iw_angle)));
        [G, patient, seizure] = findgroups(fits.patient, fits.seizure);
        
        out.patient = patient; out.seizure = seizure;
        out.iw_angle_pval = nan(max(G), 1);
        
        for mm = string(F.Metrics)
            [out.(mm).pre, out.(mm).post] = deal(nan(max(G), 2));
            for gg = unique(G)'
                dat = fits(G == gg & strcmpi(fits.metric, mm), :);
                pre = dat.dir(dat.phase == 1); pre(isnan(pre)) = [];
                post = dat.dir(dat.phase == 3); post(isnan(post)) = [];

                out.(mm).pre(gg, :) = [circ_mode(pre) numel(pre)];
                out.(mm).post(gg, :) = [circ_mode(post) numel(post)];
                out.iw_angle_pval(gg) = dat.angle_p_val(1);
            end
        end
        
    end
    
    
    function [out, mtc] = kl_divergence(F)
        % This computes a bunch of measures comparing the methods and the
        % resulting distributions; not just kl_divergence. 
        % Used to generate scatter plots
        
        RES = 64;
        
        MTC = F.Metrics;
        data = F.CompiledData;
        mask = ismember(data.metric, MTC);
        data = data(mask, :);
        
        [G, pat, sz] = findgroups(data.patient, data.seizure);
        
%         kl_dist = @(p, q) mean([ ...
%             nansum(p .* log(p./(q + eps))), ...
%             nansum(q .* log(q./(p + eps))) ...
%             ]);
        xi = linspace(-pi, pi, RES + 1);
        
        kl_dist = @(p, q) nansum(p .* log(p./q));
        
        nP = nchoosek(numel(MTC), 2);
        mpairs = nchoosek(MTC, 2);
        [Dkl, Dmode, Dmean_diffs, Dstd_diffs, Dskew_diffs, Dcorr_pval, Dcorr_coeff] = ...
            deal(nan(max(G), nP));
        for ii = 1:max(G)  % for each group
            dat = data(G == ii, :);
            [Gm, mtc] = findgroups(dat.metric);
            assert(all(ismember(MTC, mtc)));  % throw an error if something is missing
            
            
            % Compute KL-divergence based on distribution of differences.
            % Compare dist of diffs to VM normal with same center and std
            W = WaveProp.load(sprintf('%s_Seizure%d', pat{ii}, sz(ii)), F.Metrics);
            for jj = 1:nP  
                fit1 = W.(mpairs{jj, 1});
                fit2 = W.(mpairs{jj, 2});
                if numel(fit1.time) > numel(fit2.time)
                    diffs = fit1.diff(fit2);
                else
                    diffs = fit2.diff(fit1);
                end
                diffs(isnan(diffs)) = [];
                
                stats = circ_stats(diffs);
                kappa = circ_kappa(diffs);
                vm_pdf = circ_vmpdf(xi, stats.mean, kappa);
                pdfD = ksdensity([diffs - 2*pi; diffs; diffs+2*pi], xi, ...
                'bandwidth', pi/8)';
                Dkl(ii, jj) = kl_dist(pdfD./sum(pdfD), vm_pdf./sum(vm_pdf));
                Dmean_diffs(ii, jj) = stats.mean;
                Dstd_diffs(ii, jj) = stats.std;
                Dskew_diffs(ii, jj) = stats.skewness;
                
            end
            pdf = splitapply(@(x) ...
                ksdensity([x - 2*pi; x; x+2*pi], xi, ...
                'bandwidth', pi/8), dat.dir, Gm);
            
            [~, loc] = max(pdf, [], 2);
            Dmode(ii, :) = angle(exp(1j*pdist(loc, 'euclidean') * 2*pi/RES));
%             Dkl(ii, :) = pdist(pdf, kl_dist);
%             Dkl.(sprintf('%s_%d', pat{ii}, sz(ii))) = pdist(pdf, kl_dist);
            
            pairs = nchoosek(1:numel(mtc), 2);
            
            for jj = 1:nP
                p = pdf(pairs(jj, 1), :);
                q = pdf(pairs(jj, 2), :);
                
                [cc, pval] = circ_corrcc(p, q);
                Dcorr_coeff(ii, jj) = cc;
                Dcorr_pval(ii, jj) = pval;
            end
            
        end
        
        out = table(Dkl, Dmode, Dmean_diffs, Dstd_diffs, Dskew_diffs, ...
            Dcorr_pval, Dcorr_coeff, pat, sz);
        
    end
    
    function out = get.FitDiffs(F)
        if isempty(F.FitDiffs)
            W = warning; warning off;
            F.FitDiffs = F.kl_divergence;
            F.FitDiffs.Dmode_d = F.FitDiffs.Dmode / pi * 180;
            warning(W);
        end
        out = F.FitDiffs;
    end
    
    function seizure_numbers = make_polarhist(F, seizure_numbers, window)
        if nargin < 2 || isempty(seizure_numbers)
            seizure_numbers = 1:height(F.SeizureInfo);
        end
        if nargin < 3; window = []; end
        S = warning;
        warning('off');
        for ii = seizure_numbers
            [iw_info, mw] = BVNY.get_iw_info(ii);
            if isempty(iw_info), continue; end
            iw_info = iw_info{mw};
            W = WaveProp.load(F.SeizureInfo.name(ii));

            for period = ["early", "late"]
                h = F.polarhist(W, iw_info, period, window); 
                F.print(h, sprintf('%sfig3/polarhist/%s_%s', ...
                    F.Prefix, W.Name{:}, period)); 
                close(h); 
            end

        end
        warning(S);
    end

    function h = polarhist(F, W, iw_info, period, window)
        % Make a polar histogram of TW directions by seizure phase
        if nargin < 4 || isempty(period), period = ''; end
        if nargin < 5 || isempty(window), window = inf; end

        period = validatestring(period, {'Early', 'Late', '', 'All'});
        KUIPER = false;  % I ran these and it just showed that almost none of the distributions are identical
        RES = 64;

        iw_center = iw_info.onset;
        iw_t0 = iw_info.range(1);
        iw_tf = iw_info.range(2);
        iw_angle = iw_info.wave_fit.phi;

        % Create the figure
        h = figure('position', [0 0 2 2]);
        pax = polaraxes(h);
%             h.Resize = 'off';

        MTC = sort(F.Metrics);  % This sorts alphabetically (quick way to put D methods first so they are on the bottom)
        % Add data for each metric
        rmax = zeros(numel(MTC), 1);
        for ii = 1:numel(MTC)

            metric = MTC{ii};
            mtc = W.(metric);
            mtc.RotateBy = iw_angle;
            data = mtc.Direction;


            maskE = mtc.time > iw_t0 - window & mtc.time < iw_center;
            maskL = mtc.time > iw_center & mtc.time < iw_tf + window;
            non_nan = ~isnan(data);

            if KUIPER  % Use Kuiper tests to compare pre/post; metrics are compared later

                [p, k, K] = circ_kuipertest( ...
                    data(maskE & non_nan), ...
                    data(maskL & non_nan));
                [p_ww, tbl] = circ_wwtest( ...  % this may be inappropriate since it assumes von-mises and the same concentration parameter...
                    data(maskE & non_nan), ...
                    data(maskL & non_nan));
                fid = fopen(sprintf('kuiper_tests/%s_%s_EvL.txt', W.Name{:}, metric), 'w');
                fprintf(fid, 'p: %0.4g\nk: %0.4f\nK: %0.4f\np_ww: %0.4g', ...
                    p, k, K, p_ww);
                fclose(fid);

                C.(metric) = data(non_nan);  % for kuiper test on different methods
            end
            % Mask data to indicated period
            switch period

                case 'Early'
                    data = data(maskE & non_nan);
                case 'Late'
                    data = data(maskL & non_nan);
                otherwise
                    period = 'All';
            end
            if all(isnan(data)), continue; end  % CUCX5_3 has no early directions in M method

            % ... put it in a histogram
            ph = polarhistogram(pax, data, linspace(-pi, pi, RES + 1), ...
                'normalization', 'probability', ...
                'facecolor', F.Style.(metric).color, ...
                'edgecolor', F.Style.(metric).color);

            % ... prettify the plot
            set(pax, 'thetaaxisunits', 'degrees', ...
                'thetaticklabel', [], ...
                'nextplot', 'add')

            rmax(ii) = max(ph.Values);




        end

        % kuiper test on metric pairs
        if KUIPER
            pairs = nchoosek(1:numel(F.Metrics), 2);
            for ii = 1:size(pairs, 1)
                [mtc1, mtc2] = F.Metrics{pairs(ii, :)};
                [p, k, K] = circ_kuipertest(C.(mtc1), C.(mtc2));
                [p_ww, tbl] = circ_wwtest(C.(mtc1), C.(mtc2));
                fid = fopen(sprintf('kuiper_tests/%s_%sv%s.txt', W.Name{:}, mtc1, mtc2), 'w');
                fprintf(fid, 'p: %0.4g\nk: %0.4f\nK: %0.4f\np_ww: %0.4g', ...
                    p, k, K, p_ww);
                fclose(fid);
            end
        end

        rmax = ceil(rmax * 10)/10;
        pax.RDir = 'reverse';
        pax.RTick = unique([0; rmax]);
        pax.RLim = [0 max([rmax; pax.RLim(2)])];
        pax.RTickLabel{1} = '';

%             rmax = pax.RTick(end);
%             if mod(rmax, .2), rmax = rmax + .1; end
%             pax.RTick = [0 rmax/2 rmax];
%             pax.RLim(2) = max(pax.RLim(2), rmax);

        ttl = ...
            strrep( ...
                strrep( ...
                    strrep(W.Name{:}, '_fits', ''), ...
                '_', ' '), ...
            'Seizure', 'Sz.');
%             title(pax, {ttl; period})
        % hold(ax, 'off');
        pax.Tag = [ttl ' ' period];
    end
    
    function h = pre_post_direction_hist(F, interval, min_electrodes)
        % h = pre_post_direction_hist(F, interval, min_electrodes)
        
        if nargin < 2, interval = []; end
        if nargin < 3, min_electrodes = []; end
        
        out = F.transient_iw(interval, min_electrodes);
        cols =@(mm) F.Style.(mm).color;
        
        for mm = string(F.Metrics)
            out.(mm).diff = circ_dist(out.(mm).post(:, 1), out.(mm).pre(:, 1));
        end

        h = figure('units', 'inches', 'position', [0 0 3 1]*2, ...
            'name', 'pre_post_hist');
        phases = ["pre" "post" "diff"];
        for ii = 1:3
            phase = phases(ii);
            pax = polaraxes;
            subplot(1, 3 ,ii, pax)
            for mm = string(F.Metrics)
                polarhistogram(pax, out.(mm).(phase)(:, 1), linspace(-pi, pi, 13), ...
                    'facecolor', cols(mm), 'edgecolor', cols(mm));
                hold on
            end
        %     set(pax, 'rlim', [0 1]);
%             pax.RTick = [];

        %     pax.ThetaTick = -180:90:180;
            pax.ThetaAxisUnits = 'deg';
            lbl = pax.ThetaTickLabel;
            lbl(~ismember(lbl, {'0', '180', '90', '270'})) = {''};
            lbl = strrep(lbl, '270', '-90');
            lbl = strrep(lbl, '0', '0\circ');
            pax.ThetaTickLabel = lbl;
            pax.LineWidth = 1;
        %     pax.RTickLabel(1:end-1) = {''};
            title(phase)
        end

        F.print(h, F.prefix_better(''));
    end
    
    function h = pre_post_direction_stick_plots(F, interval, min_electrodes)
        % h = pre_post_direction_stick_plots(F, interval=[], min_electrodes=[])
        
        if nargin < 2, interval = []; end
        if nargin < 3, min_electrodes = []; end
        
        out = F.transient_iw(interval, min_electrodes);
        fig = @() figure('units', 'inches', 'position', [0 0 8 2.5]*.67, ...
            'name', 'pre_post_stick');
        
        h = fig();
        ax = axes(h);  

        % Note: These are the current colors assigned to M and D10. If you changes these,
        % you will need to update this
        cols = @(mm) F.Style.(mm).color;
        mask = out.iw_angle_pval < .05;

        sz = SeizureInfo;
        [~, p, pa] = findgroups(sz.patient, sz.patientAlt);
        MAP = containers.Map(p, pa);
        G = cellfun(@(x) MAP(x), out.patient);
        [~, so] = unique(G);
        
        
        Gx = [G; max(G) + 1];
        for mm = string(F.Metrics)

            xx = [Gx Gx]' + .25*[-1; 1];  % make x-values out of patient number and seizure index


            for offset = [-360 0 360]
                yy = unwrap([out.(mm).pre(:, 1) out.(mm).post(:, 1)]') / pi * 180;  % unwrap differences and convert to degrees
                yy = yy + offset;
                N = [out.(mm).pre(:, 2) out.(mm).post(:, 2)]';
                yy(N == 0) = nan;
                yy = [yy [0; 180]]; %#ok<AGROW>


                % Indicate IW and anti-IW direction
                ls = ':-:';
                theta = [-180 0 180];
                for ii = 1:numel(theta)
                    yline(theta(ii), ['k' ls(ii)], 'linewidth', 1);
                end
                hold on
                plot(ax, xx, yy, 'color', cols(mm), 'displayname', mm)
                dn = ["Pre", "Post"];
                for ii = 1:2
                plot(xx(ii, :), yy(ii, :), '.', ...
                    'color', cols(mm), 'displayname', dn(ii), ...
                    'markersize', 15);
                end


                % Indicate good IW direction fits
                gray_ = 0.1 * [1 1 1];
                plot(xx(1, mask), yy(1, mask), '.', 'color', gray_);
                plot(xx(2, mask), yy(2, mask), '.', 'color', gray_);

            end
        %     hold(ax, 'off');


        %     title(mm)
        end
        xticks(ax, unique(Gx))
        xticklabels(ax, [out.patient(so); 'HYP'])
        xtickangle(ax, 45)
        yticks(ax, -360:90:360)
        yticklabels(ax, num2str([0 90 -180 -90 0 90 180 -90 0]'))
        ylim([-1 1] * 200)
        grid on

        ln = findobj(ax, 'type', 'line');
        legend(ln([round(numel(ln)/2), end]), ... % string(F.Metrics), ...
            'location', 'eastoutside', 'box', 'off')

        F.print(h, F.prefix_better(''))

    end
    
    function h = allP_diffs_mean_v_std(F)
        
        out = F.FitDiffs;
        h = figure('units', 'inches', 'position', [0 0 3 3]);
        
        xx = out.Dmean_diffs/pi * 180;
        yy = out.Dstd_diffs/pi * 180;
        F.gscatter_pat(xx, yy, out.pat);

        hold on;
        agree = (abs(out.Dmode_d) < 90);  % & out.Dcorr_coeff > .5 
        sc = scatter(xx(~agree), yy(~agree));
        legend('off')
        set(sc, 'Marker', '*', 'MarkerEdgeColor', [0 0 0]);
        grid on
        xlabel('Mean (\circ)')
        ylabel('STD (\circ)');
        set(gca, 'fontsize', 11);

        F.print(h, F.prefix_better(''));
    end
    
    function disagree = get.Disagree(F)
        out = F.FitDiffs;
        agree = abs(out.Dmode_d) < 60 & abs(out.Dmean_diffs) < pi/4; % & out.Dcorr_coeff >= .5;
        fprintf('Disagree:\n');
        disagree = compose('%s %d', out.pat(~agree), out.sz(~agree));
        cellfun(@disp, disagree)
        
    end
    
    function allP_pdf_distance(F)
        
        out = F.FitDiffs;
        out.Dmean_diffs = rad2deg(out.Dmean_diffs);


        fields = {'Dmean_diffs', 'Dmode_d', 'Dcorr_coeff', 'Dkl'};
        labels = {'Per discharge difference (\circ)', 'Difference in modes (\circ)', ...
            'Circular correlation', 'KL divergence'};
        polar_field = [0 0 0 0];

        h = figure('name', 'pdf_distance', 'units', 'inches', ...
            'position', [0 0 1.5*(numel(fields)+1) 3] * .67);        
        

        Ts = tiledlayout(h, 1, numel(fields));
        for ii = 1:numel(fields)

            yy = out.(fields{ii});

            nexttile(Ts, ii)
            F.gscatter_pat(yy, out.pat);

            if polar_field(ii)
                ylim([-pi pi])
                yticks(-pi:pi/2:pi)
                yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
            end
            if contains(fields{ii}, 'std'), ylim([0 pi]); end
%             if strcmp(fields{ii}, 'Dmode_d')
%                 agree = (abs(yy) < 90);  % & out.Dcorr_coeff > .5 
%             end
%             if strcmp(fields{ii}, 'Dmean_diffs')
%                 hold on;
%                 sc = scatter(xx(~agree), yy(~agree));
%                 set(sc, 'Marker', '*', 'MarkerEdgeColor', [0 0 0]);
%                 hold off
%             end
            grid on
            box off
            

            ylabel(labels{ii})
            switch fields{ii}
                case 'Dmode_d'
                    ylim([-180 180]);
                    yticks(-180:60:180)
                    yline(60, 'r--', 'linewidth', 1);
                    yline(-60, 'r--', 'linewidth', 1);
                case 'Dmean_diffs'
                    ylim([-180 180]);
                    yticks(-180:90:180)
                    yline(45, 'r--', 'linewidth', 1);
                    yline(-45, 'r--', 'linewidth', 1);
                case 'Dcorr_coeff'
                    yline(.5, 'r--', 'linewidth', 1)
            end
        end

        lgd = findobj(gcf, 'type', 'legend');
        set(lgd(2:end), 'visible', 'off');
        lgd(1).Location = 'eastoutside';

%         set(Ts.Children, 'fontsize', 11)

        F.print(h, F.prefix_better(''));
    end

    function allP_main_iw(F, min_electrodes)
        if nargin < 2 || isempty(min_electrodes), min_electrodes = 0; end
        
        iw_all = F.iw_info_all;
        Nch = structfun(@(s) numfinite(s.template.template(s.main_wave, :)), iw_all);
        mask = Nch >= min_electrodes;
        h = figure('name', 'iw_templates'); fullwidth(true);
        T = tiledlayout(h, 'flow');
        
        fields = fieldnames(iw_all);
        for ff = string(fields(mask)')
            tpl = squeeze(iw_all.(ff).template.template(iw_all.(ff).main_wave, :, :));
            fr_tpl = squeeze(iw_all.(ff).template.firing_rate(iw_all.(ff).main_wave, :, :));
            [xx, yy] = find(isfinite(tpl));
            t_val = tpl(isfinite(tpl));
            fr_val = fr_tpl(isfinite(tpl));
            
            ax = nexttile(T);
            scatter(ax, xx, yy, min(fr_val*.75, 200), t_val, 'filled');
            axis square
            colorbar
            title(strrep(ff, '_', ' '))
            xlim([0 11])
            ylim([0 11])
            xticks([])
            yticks([])
            box on
        end
        
        F.print(h, F.prefix_better(''));
        
    end
    
    function allP_dir_v_rho(F, metric)
        % Correlation v. angle (relative to IW angle)
        
        data = F.CompiledData;
        mask = strcmpi(data.metric, metric) ...
            & data.nchannels >= 40; % ...
%             & data.rho_pval < 5e-2;
        data = data(mask, :);
        
        F.gscatter_pat(abs(angle(exp(1j*(data.dir-data.iw_angle)))), ...
            data.rho, data.patient);
    end
    
    function [h, sc] = show_dir_v_time_early_rasters(F, patient, metric)
        % Compare directions using all electrodes and only the first XX
        % electrodes from each discharge time
        
        if nargin < 2 || isempty(patient), patient = "c7"; end % c7
        if nargin < 3 || isempty(metric), metric = 'M'; end
        if isnumeric(patient), patient = string(F.SeizureInfo.patient(patient)); end
        
        % Get all seizures for the given patient
        sz_num = find(strcmpi(F.SeizureInfo.patient, patient));
        
        h = figure('name', 'dir raster original v first 20', ...
            'position', [0 0 2 numel(sz_num)]*1);
        T = tiledlayout(h, numel(sz_num), 1);
        
        for ss = sz_num(:)'
            
            % Load the original fit and the IW-only fit
            M = WaveProp.load(F.get_file(ss), metric);

            % Refit with only the first 20 times
            data = M.(metric).TOA;
            [~, so] = sort(fillmissing(data(:, :), 'constant', -99), 2, 'desc');
            for ii = 1:size(so, 1)
                data(ii, so(ii, 21:end)) = nan;
            end
            data(data == -99) = nan;
            M.Mr = copy(M.(metric));
            M.Mr.TOA = data;
            M.Mr = M.Mr.refit_data;
            M.Mr.MinFinite = 5;

            M0 = M.(metric);
            Mnew = M.Mr;
            dir_diff = Mnew.diff(M0);

            ax = nexttile(T);
            Mnew.direction_raster(ax, dir_diff);
            
            
        end
        xlabel(ax, 'Time [s]');
        
        sc = findobj(h, 'type', 'scatter');
        ax = findobj(h, 'type', 'axes');
        
        linkaxes(ax, 'x');
        
        tag = sprintf('%s_%s', patient, metric);
        F.print(F.prefix_better(tag));

    end
    
    function [h, sc] = show_dir_v_time_early(F, patient, metrics, style)
        % Compare directions using all electrodes and only the first XX
        % electrodes from each discharge time
        
        if nargin < 2 || isempty(patient), patient = "c7"; end % c7
        if nargin < 3 || isempty(metrics), metrics = F.Metrics; end
        if nargin < 4 || isempty(style), style = 'raster'; end
        if isnumeric(patient), patient = string(F.SeizureInfo.patient(patient)); end
        style = validatestring(style, ["raster" "scatter"]);
        
        % Get all seizures for the given patient
        sz_num = find(strcmpi(F.SeizureInfo.patient, patient));
        
        h = figure('name', 'dir scatter original v first 20', ...
            'position', [0 0 2*numel(metrics) numel(sz_num)]*1, ...
            'resize', 'off');
        T = tiledlayout(h, numel(sz_num), numel(metrics), 'padding', 'none');
        
        cmap = make_diverging_colormap([1 0 0; 1 0 0], [1 1 1], 256);
        color = @(mm) F.Style.(mm).color;
        for ss = sz_num(:)'
            % Get the IW info
            [iw, mw] = BVNY.get_iw_info(F.get_fname(ss));
            if isempty(iw), iw = struct; iw.phi = 0; else, iw = iw{mw}; end

            % Load the original fit and the IW-only fit
            M = WaveProp.load(F.get_file(ss), metrics);

            for mm = metrics
                % Refit with only the first 20 times
                toa20 = M.(mm).toa_from_first_N_electrodes(20);
                
                M.Mr = copy(M.(mm));
                M.Mr.TOA = toa20;
                M.Mr = M.Mr.refit_data;
                M.Mr.MinFinite = 5;
                
                % Rotate by the iw angle
                for ff = [mm, "Mr"], M.(ff).RotateBy = iw.phi; end
                M0 = M.(mm);
                Mnew = M.Mr;

                ax = nexttile(T);
                switch style
                    case 'scatter'
                        % Make the scatters
                        sc1 = M0.direction_scatter(ax, 'markerfacecolor', [0 0 0], 'sizedata', 2);
                        hold(ax, 'on');
                        sc2 = Mnew.direction_scatter(ax, 'markerfacecolor', [1 0 0], 'sizedata', 2);
                        hold(ax, 'off');
                        iw_dot = findobj(ax, 'tag', 'iw_dot');
                        set(iw_dot, 'markerfacecolor', color(mm), 'sizedata', 36);
                        ylabel(ax, 'TW dir');
                        set(ax.Title, 'string', strrep(ax.Title.String, 'Seizure', 'Sz.'));
                        
                    case 'raster'
                        % Make the rasters
                        dir_diff = Mnew.diff(M0);
                        Mnew.direction_raster(ax, dir_diff, cmap);
                        ylabel(ax, ["\phi_{early} - \phi_{all}"])
                end
                ax.Tag = mm;
            end
        end
        
        % Prettify
        ax = findobj(T, 'type', 'axes');
        grid(ax, 'on');
        
        % xlabel the bottom axes
        bottom_ax = ax(1:numel(metrics));
        arrayfun(@(ax) xlabel(ax, 'Time [s]'), bottom_ax); 
        % move titles to ylabels (commented out for now but axes look so
        % busy... so leave this here in case you want to change the
        % labeling and it's useful)
%         ttl2ylabel =@(ax) ...
%             {set(ax.YLabel, 'String', strsplit(strrep(ax.Title.String, 'Seizure', 'Sz.')), ...
%                 'rotation', 0, 'horizontalalignment', 'right', 'verticalalignment', 'middle'); ...
%             set(ax.Title, 'String', '')};
%         arrayfun(@(ax) ttl2ylabel(ax), ax, 'uni', 0);
        % only show ylabels on left axes
        arrayfun(@(ax) set(ax.YLabel, 'visible', 'off'), ax);
        left_ax = findobj(ax, 'tag', metrics(1));
        arrayfun(@(ax) set(ax.YLabel, 'visible', 'on'), left_ax);
        
        
        % Gather scatter objects for return
        sc = findobj(h, 'type', 'scatter');
        
        linkaxes(ax, 'x');
        
        tag = sprintf('%s_%s_%s_%s', style, patient, metrics);
        F.print(h, F.prefix_better(tag));
        
        if style == "raster"
            pax = WaveProp.colorwheel(cmap);
            F.print(pax.Parent, F.prefix_better('colorwheel'));
            close(pax.Parent);
        end
    end
    
    function sc = show_dir_v_time_iw_electrodes(F, sz_num, metric)
        % Compare directions using all electrodes and IW-only electrodes
        
        sz_num = 9;  % c7
        metric = 'M';
        metrics = {metric, [metric 'iw']};
        
        
        h = figure('name', 'dir scatter original v iw', 'position', [0 0 2*numel(sz_num) numel(sz_num)]*1);
        T = tiledlayout(h, 1, numel(sz_num));
        
        % Get the IW info
        [iw, mw] = BVNY.get_iw_info(F.get_fname(sz_num));
        if isempty(iw), iw = struct; iw.phi = 0; else, iw = iw{mw}; end
        
        % Load the original fit and the IW-only fit
        M = WaveProp.load(F.get_file(sz_num), metrics);
        for ff = string(metrics), M.(ff).RotateBy = iw.phi; end
        M0 = M.M;
        Mnew = M.Miw;
        
        ax = nexttile(T);
        sc1 = M0.direction_scatter(ax, 'markerfacecolor', [0 0 0], 'sizedata', 2);
        hold(ax, 'on');
        sc2 = Mnew.direction_scatter(ax, 'markerfacecolor', [1 0 0], 'sizedata', 2);
        hold(ax, 'off');
        
        sc = get(ax, 'children');
        
    end
    
    function sc = allP_corr_v_time(F, metric, thresh, min_electrodes)
        % Shows the correlation v time for each patient (different seizures
        % identified by color)
        % sc = allP_corr_v_time(F, metric='M', thresh=5e-2, min_electrodes=0)
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 40; end
        

        % Get the data
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) ...  % mask by metric,
            & isfinite(data.rho) ...  % ... finite rho (corr coeff; if this isn't finite, it's probably because the number of electrodes in the discharge was < F.MinFinite)
            & data.nchannels >= min_electrodes;  % ... and # electrodes recruited to IW
        data = data(mask, :);

        % Group by patient
        [G, pat] = findgroups(data.patient);
        r = numel(unique(pat));
        
        % Create the figure and tiled layout
        h = figure('name', 'corr_v_time', 'position', [0 0 3 r]);
        T = tiledlayout(h, r, 1);
        gray_ = .5*[1 1 1];

        % For each patient create a scatter plot of the rho value for each
        % discharge time. Show all rho in gray; highlight significant rho
        % in color; use different colors for each seizure
        for ii = 1:numel(pat)
            
            % Show all discharges in gray
            ax = nexttile(T, ii);
            xx = data.rho_time_offset(G == ii);
            yy = data.rho(G == ii);
            sc = scatter(ax, xx, yy, 6, 'filled', 'markerfacecolor', gray_);

            % highlight significant rho in color
            hold on
            cc = data.seizure(G == ii);  % color by patient
            mask = data.rho_pval(G == ii) < thresh;
            sc = gscatter(ax, xx(mask), yy(mask), cc(mask), [], [], 6, 'off');
            
            % prettify
            ylim([-1 1]);
            xlim([-inf inf])
            title(sprintf('%s', pat(ii)))
            xline(0)
            hold off
            grid on
            ylabel('\rho');
            
        end
        
        % prettify
        xlabel('Time [s]')
        linkaxes(T.Children, 'x');
        F.print(h, F.prefix_better(metric));

    end
    
    function dataR = allP_hist2d_tw_v_iw_effectsize(F, metric, thresh, min_electrodes, rho_or_dir)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a resampling procedure to estimate the effect size of the
        % pos/neg discharge rate in 2*HALFWIN second windows
        if nargin < 2 || isempty(metric), metric = "M"; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 40; end
        if nargin < 5 || isempty(rho_or_dir), rho_or_dir = 'rho'; end
        rho_or_dir = validatestring(rho_or_dir, ["rho" "dir"]);
        metric = string(metric);

        Nouter = 1e4;
        HALFWIN = 2;
        ptest = .05; % Show CI to this level
        if rho_or_dir == "rho"
            dnames = ["Corr" "Anticorr"];
        else
            dnames = ["d\phi < 90" "d\phi > 90"];
        end
        
        h = figure('name', 'hist2d_dir', 'position', [0 0 2 1] * 2);
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) ...
            & data.nchannels > min_electrodes;
        data = data(mask, :);

        trange = quantile(data.rho_time_offset, [0 1]);
        tbins = linspace(trange(1), trange(2), diff(trange)/2 + 1);  % Divide times into ~2s intervals


        % Group data by patient; get indices;
        % determine the lowest number of discharges in each
        % patient. Use this as the resample number
        G = findgroups(data.patient);
        N = min(splitapply(@numel, G, G));
        
        % Uncomment here to use direction instead of correlation
        if rho_or_dir == "rho"
            phi = sign(data.rho) .* double(data.rho_pval < thresh);
            sample_times = data.rho_time_offset;
        else
            dir = angle(exp(1j*data.dir - data.iw_angle));
            phi = double(abs(dir) <= pi/2) - double(abs(dir) > pi/2);
            sample_times = data.time - data.iw_center;
        end
        
        
        % Create a function to add noise (I set this arbitrarily...) 
        noise_for = @(inds) 1 * randn(size(inds)); 

        % Get the relevant indices of data for each patient (i.e. those of
        % all times and those of sig times)
        Iall = arrayfun(@(ii) find( (G == ii) ), 1:max(G), 'uni', 0);  % indices of all times for each patient
            
        % Initialize matrices to store the results
        dataR = nan(numel(tbins), 3, Nouter);
        for ii = 1:Nouter

            % Get resampled indices (bootstrap Ninner times from each patient)
            subs = cellfun(@(x) randi(numel(x), 1, N), Iall, 'uni', 0);  % random resampling of indices
            IR = cell2mat(cellfun(@(x, s) x(s), Iall, subs, 'uni', 0));  % resampled indices

            dirR = phi(IR(:));
            timeR = sample_times(IR(:)) + noise_for(IR(:));

            dir_pos = arrayfun(@(tt) sum(abs(tt - timeR(dirR == 1)) < HALFWIN), tbins);
            dir_neg = arrayfun(@(tt) sum(abs(tt - timeR(dirR == -1)) < HALFWIN), tbins);
            dir_all = arrayfun(@(tt) sum(abs(tt - timeR) < HALFWIN), tbins);
            dataR(:, :, ii) = [dir_pos; dir_neg; dir_all]';
        end
        proportions = dataR(:, [1 2], :) ./ dataR(:, 3, :);

        
        % Find time points where there are too many nans ( more than ptest)
        pct_finite = mean(isfinite(proportions), 3);
        mask = pct_finite < 1-ptest;
        
        Q = quantile(proportions, [ptest/2 1-ptest/2], 3);
        Qlo = Q(:, :, 1); Qlo(mask) = nan;
        Qhi = Q(:, :, 2); Qhi(mask) = nan;

        md = nanmedian(proportions, 3); md(mask) = nan;
        mn = nanmean(proportions, 3); mn(mask) = nan;
        sd = nanstd(proportions, [], 3); sd(mask) = nan;
        

        % *** Show the result ***
        T = tiledlayout(h, 2, 1);
        interp_ = @(x) interp(x, 5);  % interp for pretty
        fill_ = @(x) fillmissing(x, 'linear', 'endvalues', 'none');
        fillQ_ = @(x, val) fillmissing(x, 'constant', val, 'endvalues', 'none');
        
        tbinsI = interp_(tbins);
        
        
        % Show median and 95CI
        mnI = [interp_(fill_(md(:, 1))) interp_(fill_(md(:, 2)))];
        QloI = [interp_(fillQ_(Qlo(:, 1), 0)) interp_(fillQ_(Qlo(:, 2), 0))];
        QhiI = [interp_(fillQ_(Qhi(:, 1), 1)) interp_(fillQ_(Qhi(:, 2), 1))];
        
        tt = [tbinsI fliplr(tbinsI)]';

        
        m0 = strsplit(metric, 'sub');
        color = F.Style.(m0(1)).color;
        gray_ = .5 * [1 1 1];
        for jj = 1:2

            % Show the CI
            ax = nexttile(T, jj);
            
%             yy = fill_(mnI(:, jj) + 3*[-1 1] .* sdI(:, jj));  % sd limits
            yy = [QloI(:, jj) QhiI(:, jj)];  % quantile limits
            
            yyF = [yy(:, 1); flipud(yy(:, 2))];
            mask = isfinite(yyF);
            temp = fill(ax, tt(mask), yyF(mask), 1, ...
                'facecolor', color, 'facealpha', .5, 'linestyle', 'none'); %#ok<NASGU>

            % highlight where lower bound is above 0
            hold(ax, 'on');
            mask = yy(:, 1) > 0;
            temp = area(ax, tbinsI(mask), yy(mask, 1), ...
                'facecolor', gray_, 'linestyle', 'none'); %#ok<NASGU>

            % Show the mean
            plot(ax, tbinsI, mnI(:, jj), 'color', color, ...
                'linewidth', 2, 'displayname', dnames(jj))
            hold(ax, 'off');

            % highlight IW time 
            yline(ax, 0)

            % Prettify
            grid(ax, 'on')
            axis(ax, 'tight');
            title(ax, dnames(jj));

        end

        linkaxes(T.Children, 'xy');
        ylim(ax, [0 1]);
%         xlim(ax, [-30 50]);  % Nothing significant beyond this
        xlabel(ax, 'Time [s]');
        ylabel(T, 'Rate');
        temp = legend({'95%CI', 'CI>0', 'median'}, 'location', 'eastoutside'); %#ok<NASGU>


        % print result
        F.print(h, F.prefix_better(sprintf('%s_%s', metric, rho_or_dir)));
    end
    
    
    function dataR = ZZallP_hist2d_rho_tw_v_iw_effectsize(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a resampling procedure to estimate the effect size of the
        % pos/neg discharge rate in 2*HALFWIN second windows
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 40; end

        Nouter = 1e4;
        HALFWIN = 2;

        h = figure('name', 'hist2d_rho', 'position', [0 0 2 1] * 2);
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) ...
            & data.nchannels > min_electrodes;
        data = data(mask, :);

        trange = quantile(data.rho_time_offset, [0 1]);
        tbins = linspace(trange(1), trange(2), diff(trange)/2 + 1);  % Divide times into ~2s intervals


        % Group data by patient; get indices;
        % determine the lowest number of discharges in each
        % patient. Use this as the resample number
        G = findgroups(data.patient);
        N = min(splitapply(@numel, G, G));
        rho = sign(data.rho) .* double(data.rho_pval < thresh);
        
        % Uncomment here to use direction instead of correlation
%         dir = angle(exp(1j*data.dir - data.iw_angle));
%         rho = double(abs(dir) <= pi/2) - double(abs(dir) >= pi/2);

        
        % Create a function to add noise (I set this arbitrarily...) and
        % function that will shuffle a set of indices
        noise_for = @(inds) 1 * randn(size(inds)); 

        % Get the relevant indices of data for each patient (i.e. those of
        % all times and those of sig times)
        Iall = arrayfun(@(ii) find( (G == ii) ), 1:max(G), 'uni', 0);  % indices of all times for each patient

        
        % Get our observation
        obs_pos = arrayfun(@(tt) sum(abs(tt - data.rho_time_offset(rho == 1)) < HALFWIN), tbins);
        obs_neg = arrayfun(@(tt) sum(abs(tt - data.rho_time_offset(rho == -1)) < HALFWIN), tbins);
        obs_all = arrayfun(@(tt) sum(abs(tt - data.rho_time_offset) < HALFWIN), tbins);
        
        obs = ([obs_pos; obs_neg] ./ obs_all)';
            
        % Initialize matrices to store the results
        dataR = nan(numel(tbins), 3, Nouter);
        for ii = 1:Nouter

            % Get resampled indices (bootstrap Ninner times from each patient)
            subs = cellfun(@(x) randi(numel(x), 1, N), Iall, 'uni', 0);  % random resampling of indices
            IR = cell2mat(cellfun(@(x, s) x(s), Iall, subs, 'uni', 0));  % resampled indices

            rhoR = rho(IR(:));
            timeR = data.rho_time_offset(IR(:)) + noise_for(IR(:));

            rho_pos = arrayfun(@(tt) sum(abs(tt - timeR(rhoR == 1)) < HALFWIN), tbins);
            rho_neg = arrayfun(@(tt) sum(abs(tt - timeR(rhoR == -1)) < HALFWIN), tbins);
            rho_all = arrayfun(@(tt) sum(abs(tt - timeR) < HALFWIN), tbins);
            dataR(:, :, ii) = [rho_pos; rho_neg; rho_all]';
        end
        ybins = [1 -1];
        proportions = dataR(:, [1 2], :) ./ dataR(:, 3, :);

        ptest = .05; % move this up
        
        % Find time points where there are too many nans ( more than ptest)
        pct_finite = mean(isfinite(proportions), 3);
        mask = pct_finite < 1-ptest;
        
        Q = quantile(proportions, [ptest/2 1-ptest/2], 3);
        Qlo = Q(:, :, 1); Qlo(mask) = nan;
        Qhi = Q(:, :, 2); Qhi(mask) = nan;
%         Q = min(proportions, [], 3);
        md = nanmedian(proportions, 3); md(mask) = nan;
        mn = nanmean(proportions, 3); mn(mask) = nan;
        sd = nanstd(proportions, [], 3); sd(mask) = nan;
        z = mn ./ sd;


        % *** Show the result ***
        interp_ = @(x) interp(x, 5);
        fill_ = @(x) fillmissing(x, 'linear', 'endvalues', 'none');
        fillQ_ = @(x, val) fillmissing(x, 'constant', val, 'endvalues', 'none');
        
        tbinsI = interp_(tbins);
        
        % Show mean and std
        mnI = [interp_(fill_(mn(:, 1))) interp_(fill_(mn(:, 2)))];
        sdI = [interp_(fill_(sd(:, 1))) interp_(fill_(sd(:, 2)))];
        
        % Show median and 95CI
        mnI = [interp_(fill_(md(:, 1))) interp_(fill_(md(:, 2)))];
        QloI = [interp_(fillQ_(Qlo(:, 1), 0)) interp_(fillQ_(Qlo(:, 2), 0))];
        QhiI = [interp_(fillQ_(Qhi(:, 1), 1)) interp_(fillQ_(Qhi(:, 2), 1))];
        
        T = tiledlayout(h, 2, 1);
        
        % show normalized difference
        tt = [tbinsI fliplr(tbinsI)]';
%         yy_pos = [Q(:, 1, 1); flipud(Q(:, 1, 2))];  % quantiles instead of std

        dnames = ["Corr" "Anticorr"];
        m0 = strsplit(metric, 'sub');
        color = F.Style.(m0(1)).color;
        gray_ = .5 * [1 1 1];
        for jj = 1:2

            % Show the CI
            ax = nexttile(T, jj);
            
%             yy = fill_(mnI(:, jj) + 3*[-1 1] .* sdI(:, jj));  % sd limits
            yy = [QloI(:, jj) QhiI(:, jj)];  % quantile limits
            
            yyF = [yy(:, 1); flipud(yy(:, 2))];
            mask = isfinite(yyF);
            temp = fill(ax, tt(mask), yyF(mask), 1, ...
                'facecolor', color, 'facealpha', .5, 'linestyle', 'none');

            % highlight where lower bound is above 0
            hold(ax, 'on');
            mask = yy(:, 1) > 0;
            temp = area(ax, tbinsI(mask), yy(mask, 1), ...
                'facecolor', gray_, 'linestyle', 'none');

            % Show the mean
            plot(ax, tbinsI, mnI(:, jj), 'color', color, ...
                'linewidth', 2, 'displayname', dnames(jj))
            hold(ax, 'off');

            % highlight IW time and effect size = 0
%             xline(ax, 0)
            yline(ax, 0)

            % Prettify
            grid(ax, 'on')
            axis(ax, 'tight');
            title(ax, dnames(jj));

        end

        linkaxes(T.Children, 'xy');
        ylim(ax, [0 1]);
%         xlim(ax, [-30 50]);  % Nothing significant beyond this
        xlabel(ax, 'Time [s]');
        ylabel(T, 'Rate');
        temp = legend({'95%CI', 'CI>0', 'median'}, 'location', 'eastoutside');



        % imagesc visualization
%         imagesc(tbins, ybins, z', [3 inf]); axis xy
%         ax = gca;
%         grid(ax, 'on')
%         colormap(ax, 1-gray)
%         cb = colorbar;
%         title(cb, ["Z" "(\mu/\sigma)"]);
% 
%         % show all points where difference is significant to ptest threshold
%         [xmask, ymask] = find(Q > 0);      
%         hold on; plot(tbins(xmask), .5*ybins(ymask), 'r.'); hold off
% 
%         % highlight IW crossing time (t=0)
%         xline(0)
%         yline(0)
% 
%         % prettify
%         legend(sprintf('\\lambda(t)>0, p<%4.2g', 1/Nouter), 'box', 'off', 'location', 'northoutside')
%         title(sprintf('Effect size \\lambda(t) (%s)', metric))
%         xlabel('Time [s]');
%         ylabel(sprintf('\\rho (p<5e%0.0f)', log10(thresh/5)))
%     %     xlim([-30 30]) 
%         yticks([-1 0 1])
%         ylim([-1 1])


        % print result
        F.print(h, F.prefix_better(metric));
    end

    function [counts, shuffles] = allP_hist2d_rho_tw_v_iw_v0(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a bootstrapping & shuffling procedure. 
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
        Nouter = 1e4;
        Ninner = 100;
        
        h = figure('name', 'hist2d_rho', 'position', [0 0 1.5 1] * 2);
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) & data.nchannels > min_electrodes;
        data = data(mask, :);

        trange = quantile(data.rho_time_offset, [0 1]);
        xbins = linspace(trange(1), trange(2), diff(trange)/2 + 1);  % Divide times into ~2s intervals
        ybins = linspace(-1, 1, 3);  % divide correlations into positive and negative

        
        % Group data by patient; get indices where rho_pval > thresh;
        % determine the lowest number of significant discharges in each
        % patient. Use this as the bootstrap resample number
        G = findgroups(data.patient);
        sig_mask = data.rho_pval < thresh;
        N = min(splitapply(@numel, G(sig_mask), G(sig_mask)));
       
        % Create a function to add noise (I set this arbitrarily...) and
        % function that will shuffle a set of indices
        noise_for = @(inds) 1 * randn(size(cat(1, inds{:})));
        shuffle = @(X) X(randperm(numel(X)));  
        
        % Get the relevant indices of data for each patient (i.e. those of
        % all times and those of sig times)
        I = arrayfun(@(ii) find( (G == ii) & sig_mask ), 1:max(G), 'uni', 0);  % indices of high sig
        Iall = arrayfun(@(ii) find( (G == ii) ), 1:max(G), 'uni', 0);  % indices of all times for each patient

        % Initialize matrices to store the results
        [counts, shuffles] = deal(nan(numel(xbins) - 1, numel(ybins) - 1, Nouter)); 


        for ii = 1:Nouter

            % Get resampled indices (bootstrap Ninner times from each patient)
            subs = cellfun(@(x) randi(numel(x), 1, N * Ninner), I, 'uni', 0);  % random resampling of indices
            IR = cellfun(@(x, s) x(s), I, subs, 'uni', 0);  % resampled high-sig indices

            % Shuffle the times (get a random permutation of the times
            % across all patients... is this right?)
%             t_shuff = arrayfun(@(ii) shuffle(data.rho_time_offset(Iall{ii})), 1:max(G), 'uni', 0);
%             t_shuff = cat(1, t_shuff{:});
            t_shuff = shuffle(data.rho_time_offset);

            % Get resampled rho values
            yy = data.rho(cat(1, IR{:}));

            % Create pdfs using observed times
            xx = data.rho_time_offset(cat(1, IR{:})) + noise_for(IR);
            counts(:, :, ii) = histcounts2(xx, yy, xbins, ybins, 'normalization', 'pdf');

            % Create pdfs using shuffled times
            xx = t_shuff(cat(1, IR{:})) + noise_for(IR);  % get shuffled times
            shuffles(:, :, ii) = histcounts2(xx, yy, xbins, ybins, 'normalization', 'pdf');
        end

        ptest = max(1e-4, 1/Nouter);  % set the test threshold
        Q_shuff = quantile(shuffles, 1-ptest, 3);  % get upper bound of the shuffles
        Q_data = quantile(counts, ptest, 3);  % ... and lower bound of the data
        shuff_mad = mad(shuffles, [], 3);  % Compute the MAD of the shuffles

        
        
        % Get the difference between the lower bound of the data samples
        % and the upper bound of the shuffles; normalize the difference by
        % the MAD of the shuffles.
        Qdiff = Q_data - Q_shuff;  
        Qn = Qdiff ./ shuff_mad;  
        

        
        % *** Show the result ***
        
        % show normalized difference
        imagesc(xbins, ybins, Qn', [0 inf]); axis xy; 
        ax = gca;
        grid(ax, 'on')
        colormap(ax, 1-gray)
        cb = colorbar;
        cb.Label.String = 'MAD';
        
        % show all points where difference is significant to ptest threshold
        [xmask, ymask] = find(Qdiff > 0);  
        yc = (ybins(1:end-1) + ybins(2:end))/2;
        xc = movmean(xbins, 2); xc = xc(2:end);
        hold on; plot(xc(xmask), yc(ymask), 'r.'); hold off
        
        % highlight IW crossing time (t=0)
        xline(0)
        
        % prettify
        legend(sprintf('p<%4.2g', ptest), 'box', 'off', 'location', 'northoutside')
        title(sprintf('PDF \\rho (%s, p<5e%0.0f)', metric, log10(thresh/5)))
        xlabel('Time [s]');
        xlim([-30 30])  % No significant values outside of this. If you change things, you might want to double check this.
        yticks([-1 0 1])
        ylim([-1 1])
        ylabel('CC')
        
        % print result
        F.print(h, F.prefix_better(metric));
    end

    function mdl = allP_glm_rho_tw_v_iw(F, metric, thresh, min_electrodes)
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
        data = F.CompiledData;
        mask = ...
            strcmp(data.metric, metric) ...  % limit to metric
            & isfinite(data.rho) ...  % ... and correlation can be estimated
            & data.nchannels >= min_electrodes;  % ... and IW observed on at least min_electrodes
        data = data(mask, :);
        
        trange = quantile(data.rho_time_offset, [0 1]);
        tbins = linspace(trange(1), trange(2), diff(trange)/5 + 1);  % Divide times into intervals


        % Need: #(rho > 0) for each time bin; random effects: patient, seizure
        % N ~ bin
        
%         data.bin = categorical(discretize(data.rho_time_offset, tbins));
        data.rho_pos = data.rho > 0 & data.rho_pval < thresh;
        data.rho_neg = data.rho < 0 & data.rho_pval < thresh;
        [G, pat] = findgroups(data.patient);
        
        T = table();
        tbins = -20:40;
        for tt = tbins
            temp = table();
            Npos = splitapply(@sum, abs(data.rho_time_offset - tt) < 2 & data.rho_pos, G);
            Nneg = splitapply(@sum, abs(data.rho_time_offset - tt) < 2 & data.rho_neg, G);
            temp.pat = pat;
            temp.Npos = Npos;
            temp.Nneg = Nneg;
            temp.bin = tt * ones(size(Nneg));
            T = [T; temp]; %#ok<AGROW>
        end
        T.bin = categorical(T.bin);
%         T0 = T;
%         T = T0(T0.pat ~= "CUCX3", :);
        mPos = fitglme(T, 'Npos ~ time_bin + (1|patient)', 'Distribution', 'poisson'); 
        mNeg = fitglme(T, 'Nneg ~ bin + (1|pat)', 'Distribution', 'poisson');
        est = [mNeg.Coefficients.Estimate, mPos.Coefficients.Estimate];
        ci = [mNeg.coefCI, mPos.coefCI];
        lb = ci(:, [1 3]); ub = ci(:, [2 4]);
        
        xx = [tbins(:); flipud(tbins(:))];
        yy_p = [ci(:, 3); flipud(ci(:, 4))];
        yy_n = [ci(:, 1); flipud(ci(:, 2))];
        cc = ones(size(yy_p)) .* [1 5];
        P = fill([xx xx], max([yy_p yy_n], 0) .* [1 -1], cc, 'facealpha', .5, 'linestyle', 'none');
        hold on; plot(tbins, max(est, 0) .* [-1 1]); hold off
        legend('corr', 'anti-corr')
        ax = gca;
        lbl = get(ax, 'yticklabel');
        ax.YTickLabel = strrep(lbl, '-', '');
        ylabel('\lambda')
        xlabel('Time [s]')
%         imagesc(tbins, [-1 1], lb', [0 inf]); colormap(1-gray); axis xy; colorbar
        
        
    end
    
    function [obs, shuffles] = allP_hist2d_rho_tw_v_iw(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a bootstrapping & shuffling procedure. 
        % Shuffle from the rho values. Keep times fixed
        
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
Nouter = 1e4;
ptest = .01;

h = figure('name', 'hist2d_rho', 'position', [0 0 1.5 1] * 2);
data = F.CompiledData;
mask = ...
    strcmp(data.metric, metric) ...  % limit to metric
    & isfinite(data.rho) ...  % ... and correlation can be estimated
    & data.nchannels >= min_electrodes;  % ... and IW observed on at least min_electrodes
data = data(mask, :);

trange = quantile(data.rho_time_offset, [0 1]);
xbins = linspace(trange(1), trange(2), diff(trange)/2 + 1);  % Divide times into intervals

% Make rho into {-1, 0, 1} and get groups by patient
[G, pat] = findgroups(data.patient);
rho = sign(data.rho) .* double(data.rho_pval < thresh);
rho_time = data.rho_time_offset;

% Function to shuffle the data (shuffle the indices)
inds0 = arrayfun(@(ii) find(G == ii), 1:max(G), 'uni', 0);
shuffle_fun =@() cellfun(@(x) x(randperm(numel(x))), inds0, 'uni', 0);  % get a shuffle of indices
bs_fun = @(tt, mask) histcounts(tt(mask), xbins);  % get counts in each interval


% Get the observation: proportion of pos/neg rho in each bin
pos_count = sum(splitapply(@(tt, mask) bs_fun(tt, mask), rho_time, rho > 0, G));
neg_count = sum(splitapply(@(tt, mask) bs_fun(tt, mask), rho_time, rho < 0, G));
obs = [neg_count; pos_count];


% Get the H0 distribution: resample rho values in each patient;
% repeat Nouter times

% Initialize the result matrix
shuffles = nan(2, numel(xbins) - 1, Nouter);
% rho_time2 = rho_time(cat(1, inds0{:}));
for ii = 1:Nouter

    % Shuffle the rho values in each patient
    inds_shuff = shuffle_fun();
%     rhoR = cell2mat(splitapply(@(x) {shuffle_fun(x)}, rho, G));

    % Get the resampled observation: proportion of pos/neg rho in each bin
    posR = cellfun(@(i0, iSh) bs_fun(rho_time(i0), rho(iSh) > 0), inds0, inds_shuff, 'uni', 0);
    negR = cellfun(@(i0, iSh) bs_fun(rho_time(i0), rho(iSh) < 0), inds0, inds_shuff, 'uni', 0);

    posR = sum(cat(1, posR{:}));
    negR = sum(cat(1, negR{:}));
    
    shuffles(1, :, ii) = negR; 
    shuffles(2, :, ii) = posR; 
    
end

% Q = quantile(shuffles, 1-ptest, 3);
pctl = mean(obs > shuffles, 3);

tt = movmean(xbins, 2); tt = tt(2:end); % get bin centers
imagesc(tt, [-1 1], pctl, [.95 1]); axis xy

ax = gca;
grid(ax, 'on')
colormap(ax, 1-gray)
cb = colorbar;
cb.Label.String = 'Percentile';

% show all points where difference is significant to ptest threshold
inds = find(pctl > 1-ptest);  
tnew = [tt; tt];
val = ones(size(tnew)) .* [-1; 1]*.5;
hold on; plot(tnew(inds), val(inds), 'r.'); hold off


% prettify
xline(0)
legend(sprintf('p<%4.2g', ptest), 'box', 'off', 'location', 'northoutside')
title(sprintf('PDF \\rho (%s, p<5e%0.0f)', metric, log10(thresh/5)))
xlabel('Time [s]');
% xlim([-30 30])  % No significant values outside of this. If you change things, you might want to double check this.
yticks([-1 0 1])
ylim([-1 1])
ylabel('CC')

% print result
F.print(h, F.prefix_better(metric));

    end

    function [obs, shuffles] = allP_hist2d_rho_tw_v_iw_v2(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a bootstrapping & shuffling procedure. 
        % Shuffle from the rho values. Keep times fixed
        
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
Nouter = 1e4;
ptest = .01;

h = figure('name', 'hist2d_rho', 'position', [0 0 1.5 1] * 2);
data = F.CompiledData;
mask = ...
    strcmp(data.metric, metric) ...  % limit to metric
    & isfinite(data.rho) ...  % ... and correlation can be estimated
    & data.nchannels >= min_electrodes;  % ... and IW observed on at least min_electrodes
data = data(mask, :);

trange = quantile(data.rho_time_offset, [0 1]);
xbins = linspace(trange(1), trange(2), diff(trange)/2 + 1);  % Divide times into intervals

% Make rho into {-1, 0, 1} and get groups by patient
[G, pat] = findgroups(data.patient);
rho = sign(data.rho) .* double(data.rho_pval < thresh);
rho_time = data.rho_time_offset;

% Function to shuffle the data (shuffle the indices)
inds0 = arrayfun(@(ii) find(G == ii), 1:max(G), 'uni', 0);
shuffle_fun =@() cellfun(@(x) x(randperm(numel(x))), inds0, 'uni', 0);  % get a shuffle of indices
bs_fun = @(tt, mask) histcounts(tt(mask), xbins) ./ histcounts(tt, xbins);  % get counts in each interval


% Get the observation: proportion of pos/neg rho in each bin
pos_count = nanmean(splitapply(@(tt, mask) bs_fun(tt, mask), rho_time, rho > 0, G));
neg_count = nanmean(splitapply(@(tt, mask) bs_fun(tt, mask), rho_time, rho < 0, G));
obs = [neg_count; pos_count];


% Get the H0 distribution: resample rho values in each patient;
% repeat Nouter times

% Initialize the result matrix
shuffles = nan(2, numel(xbins) - 1, Nouter);
rho_time2 = rho_time(cat(1, inds0{:}));
for ii = 1:Nouter

    % Shuffle the rho values in each patient
    inds_shuff = shuffle_fun();
%     rhoR = cell2mat(splitapply(@(x) {shuffle_fun(x)}, rho, G));

    % Get the resampled observation: proportion of pos/neg rho in each bin
    posR = cellfun(@(i0, iSh) bs_fun(rho_time(i0), rho(iSh) > 0), inds0, inds_shuff, 'uni', 0);
    negR = cellfun(@(i0, iSh) bs_fun(rho_time(i0), rho(iSh) < 0), inds0, inds_shuff, 'uni', 0);

    posR = nanmean(cat(1, posR{:}));
    negR = nanmean(cat(1, negR{:}));
    
    shuffles(1, :, ii) = negR; 
    shuffles(2, :, ii) = posR; 
    
end

% Q = quantile(shuffles, 1-ptest, 3);
pctl = mean(obs > shuffles, 3);

tt = movmean(xbins, 2); tt = tt(2:end); % get bin centers
imagesc(tt, [-1 1], pctl, [.95 1]); axis xy

ax = gca;
grid(ax, 'on')
colormap(ax, 1-gray)
cb = colorbar;
cb.Label.String = 'Percentile';

% show all points where difference is significant to ptest threshold
inds = find(pctl > 1-ptest);  
tnew = [tt; tt];
val = ones(size(tnew)) .* [-1; 1]*.5;
hold on; plot(tnew(inds), val(inds), 'r.'); hold off


% prettify
xline(0)
legend(sprintf('p<%4.2g', ptest), 'box', 'off', 'location', 'northoutside')
title(sprintf('PDF \\rho (%s, p<5e%0.0f)', metric, log10(thresh/5)))
xlabel('Time [s]');
% xlim([-30 30])  % No significant values outside of this. If you change things, you might want to double check this.
yticks([-1 0 1])
ylim([-1 1])
ylabel('CC')

% print result
F.print(h, F.prefix_better(metric));

    end

    function [obs, shuffles] = allP_hist2d_rho_tw_v_iw_v1(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a bootstrapping & shuffling procedure. 
        % Shuffle from the rho values. Keep times fixed
        
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
Nouter = 1e4;

h = figure('name', 'hist2d_rho', 'position', [0 0 1.5 1] * 2);
data = F.CompiledData;
mask = ...
    strcmp(data.metric, metric) ...  % limit to metric
    & isfinite(data.rho) ...  % ... and correlation can be estimated
    & data.nchannels >= min_electrodes;  % ... and IW observed on at least min_electrodes
data = data(mask, :);

trange = quantile(data.rho_time_offset, [0 1]);
xbins = linspace(trange(1), trange(2), diff(trange) + 1);  % Divide times into ~2s intervals

% Make rho into {-1, 0, 1} and get groups by patient
[G, pat] = findgroups(data.patient);
rho = sign(data.rho) .* double(data.rho_pval < thresh);

% Function to resample
resample = @(x) x(randi(numel(x), size(x)));


% Get the observation: mean proportion of pos/neg rho in each bin
pos_count = splitapply(@(t, rho) ... % for each patient
    histcounts(t(rho == 1), xbins), ...  % get #(rho > 0) in each time bin
    data.rho_time_offset, rho, G);  % using (time_offset, rho) for groups in G (patient)
neg_count = splitapply(@(t, rho) histcounts(t(rho == -1), xbins), data.rho_time_offset, rho, G);
all_count = splitapply(@(t) histcounts(t, xbins), data.rho_time_offset, G);

pos = mean(pos_count ./ all_count);
neg = mean(neg_count ./ all_count);
obs = [pos; neg];

% Get the H0 distribution: resample rho values in each patient;
% repeat Nouter times

% Initialize the result matrix
shuffles = nan(2, numel(xbins) - 1, Nouter);
for ii = 1:Nouter

    % Resample the rho values in each patient
    rhoR = cell2mat(splitapply(@(x) {resample(x)}, rho, G));

    % Get the resampled observation: proportion of pos/neg rho in each bin
    posR = splitapply(@(t, rho) ... % for each patient
        histcounts(t(rho == 1), xbins), ...  % get #(rho > 0) in each time bin
        data.rho_time_offset, rhoR, G);  % using (time_offset, rho) for groups in G (patient)
    negR = splitapply(@(t, rho) histcounts(t(rho == -1), xbins), data.rho_time_offset, rhoR, G);


    shuffles(1, :, ii) =  mean(posR ./ all_count);
    shuffles(2, :, ii) = mean(negR ./ all_count);
end

Q = quantile(shuffles, .95, 3);
tt = movmean(xbins, 2); tt = tt(2:end);

plot(tt, max(obs - Q, 0));
legend('pos', 'neg')

    end
    
    function [counts, counts_adj] = ZZallP_hist2d_rho_tw_v_iw(F, metric, thresh, min_electrodes)
        % Shows a 2d histogram of the correlations between the tw and the
        % iw. Uses a bootstrapping procedure. 
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 0; end
        
        h = figure('name', 'hist2d_rho', 'position', [0 0 1.5 1] * 2);
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) & data.nchannels > min_electrodes;
        data = data(mask, :);


        trange = quantile(data.rho_time_offset, [0 1]);
        xbins = linspace(trange(1), trange(2), diff(trange) + 1);
        ybins = linspace(-1, 1, 3);

        % Resample N points from each patient to generate a distribution, then
        % repeat this Nouter times
        G = findgroups(data.patient);
        sig_mask = data.rho_pval < thresh;
        N = min(splitapply(@numel, G(sig_mask), G(sig_mask)));
        disp(N)
        Nouter = 1e2;
        inds = deal(nan(max(G), N));
        
        for dummy = 1:Nouter
            for ii = 1:max(G)  % Resample N points from each patient
                
                % Get high sig indices
                X = find( (G == ii) & sig_mask );
                subs = randi(numel(X), 1, N);
                inds(ii, :) = X(subs);
                
            end
            
            % Get rhos that are significantly different from 0
            xx = data.rho_time_offset(inds(:));
            yy = data.rho(inds(:));

            % compute the pdf
            if dummy == 1
                counts = histcounts2(xx, yy, xbins, ybins, 'normalization', 'pdf');
            else
                counts_new = histcounts2(xx, yy, xbins, ybins, 'normalization', 'pdf');
                counts = counts + counts_new;
            end
        end

        % take the average of the resulting pdf
        counts = counts/dummy;
        shuffles = nan([size(counts), Nouter]);
        for ii = 1:numel(xbins)*1000
            pp = randperm(size(counts, 1));
            shuffles(:, :, ii) = counts(pp, :);
        end
        
        counts_adj = counts - quantile(shuffles, .95, 3);
        
%         counts_adjSM = counts_adj;
        % Smooth over 5 seconds 
%         counts_adjSM = smoothdata(counts_adjSM, 2, 'gaussian', .2, 'samplepoints', ybins(2:end));
%         counts_adjSM = smoothdata(counts_adjSM, 1, 'gaussian', 5, 'SamplePoints', xbins(2:end));
        
        % Show the result
%         imagesc(xbins, ybins, counts_adjSM', [0 .01]); axis xy; colorbar
        imagesc(xbins, ybins, counts_adj' > 0); axis xy; 
%         [v_max, v_ind] = max(counts_adjSM, [], 2);
        
        [v_max, v_ind] = max(counts - quantile(shuffles, .99, 3), [], 2);
        lvl = max(quantile(v_max, 10/diff(trange)), 0);
        mask = v_max > lvl; 
        yc = (ybins(1:end-1) + ybins(2:end))/2;
        xc = movmean(xbins, 2); xc = xc(2:end);
        hold on; plot(xc(mask), yc(v_ind(mask)), 'r.'); hold off
        xline(0)
        ax = gca;
        grid(ax, 'on')
        colormap(ax, 1-gray)
%         legend(sprintf('Peak (pdf>%0.4g)', lvl), 'box', 'off', 'location', 'northoutside')
        legend(sprintf('p<.01'), 'box', 'off', 'location', 'northoutside')
        title(sprintf('PDF \\rho (%s, p<5e%0.0f)', metric, log10(thresh/5)))
        xlabel('Time [s]');
        yticks([-1 0 1])
        ylim([-1 1])
        ylabel('CC')
        
        F.print(h, F.prefix_better(metric));
    end
    
    function [xx, yy, pat] = allP_quantiles_pval_tw_v_iw(F, metrics, thresh, min_electrodes)
        % Plots the percent of the discharges that have strong
        % relationships with the main IW template. i.e. Percentage of
        % correlations with p-value less than thresh for each patient and
        % seizure.
        if nargin < 2 || isempty(metrics), metrics = string(F.Metrics); end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 40; end
        
        if ischar(metrics) || numel(metrics) == 1, metrics = [string(metrics) ""]; end
        
        metrics = arrayfun(@(s) validatestring(s, ["", F.Metrics]), metrics);
        figure('units', 'inches', 'position', [0 0 1.5 1] * 2, ...
            'name', 'TW-IW corr');
        
        data = F.CompiledData;
        mask = ...  % only use data where
            ismember(data.metric, metrics) ...  % the TOA method matches
            & isfinite(data.rho_pval) ...  % the correlation could be computed
            & data.nchannels >= min_electrodes;  % the IW was detected on at least min_electrodes
        data = data(mask, :);

        [G, pat] = findgroups(data.patient, data.seizure);
        
        m1 = metrics(1); m2 = metrics(2);
        mask =@(mm) strcmpi(data.metric, mm);
        if m2 == ""
%             yy = splitapply(@(x) sum(x < thresh), data.rho_pval(mask(m1)), G(mask(m1)));
            yy = splitapply(@numel, data.rho_pval(mask(m1)), G(mask(m1)));
            ly = strcat(m1, " [#]");
            tag = m1;
        else
            yy = splitapply(@(x) mean(x < thresh)*100, data.rho_pval(mask(m2)), G(mask(m2)));
            ly = strcat(m2, " [% corr\neq0]");
            tag = sprintf('%s_%s', m1, m2);
        end
        xx = splitapply(@(x) mean(x < thresh)*100, data.rho_pval(mask(m1)), G(mask(m1)));
        ln = F.gscatter_pat(xx, yy, pat);
        xlabel(strcat(m1, " [% corr\neq0]"))
        ylabel(ly)
        title({'TW detected'; ...
            sprintf('(p<5e%0.0f)', log10(thresh/5))});
        axis square
        legend('location', 'eastoutside')
%         xlim([0 100]);
%         ylim([0 100]);
        F.print(F.prefix_better(tag));
        
    end
    
    function [yy, ax] = allP_time_of_first_shift(F, metric, thresh)
        data = F.CompiledData;
        if nargin < 2 || isempty(metric), metric = 'M10'; end
        if nargin < 3, thresh = pi/2; end
        mask = strcmp(data.metric, metric) & data.dtheta_dt >= thresh;
        data = data(mask, :);
        
        [G, pat, ~] = findgroups(data.patient, data.seizure);
        yy = splitapply(@(x) x(1), data.time - data.iw_center, G);
        ax = F.gscatter_pat(yy, pat);
        title(sprintf('First shift (%s)', metric))
        ylabel('Time [s]')
        
    end
    
    function allP_proportion_finite(F, metrics)
        % This shows the number of unique seconds where each method detects
        % a TW. ::This feels a little weird so maybe there's a better way
        % but the goal is to quantify the fact that the D10 method often
        % picks up a lot of estimates that the M method doesn't. But
        % because of the differences in sampling times, it's not
        % straightforward to just compute the proportions or number of
        % detections. ::
        if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
        if numel(metrics) == 1, metrics = {'', metrics}; end
        
        
        metrics = string(metrics);
        data = F.CompiledData;
        
        figure('units', 'inches', 'position', [0 0 1.5 2]);
        mask = ismember(data.metric, metrics) & isfinite(data.dir);
        data = data(mask, :);

        m1 = metrics(1); m2 = metrics(2);

        [G, pat, sz] = findgroups(data.patient, data.seizure);

        mask =@(mm) strcmpi(data.metric, mm);

        xx = splitapply(@(x) numel(unique(round(x))), data.time(mask(m1)), G(mask(m1)));
        yy = splitapply(@(x) numel(unique(round(x))), data.time(mask(m2)), G(mask(m2)));
        F.gscatter_pat(xx, yy, pat);
        
        disagree = F.Disagree;
        dd = cellfun(@(x) find(strcmpi(x, compose('%s %d', pat, sz))), disagree);
        
        ax = gca;
        hold(ax, 'on')
        plot(ax, xx(dd), yy(dd), 'k*', 'markersize', 6);
        hold(ax, 'off')
        
        title(["Unique seconds" "[#]"]);
        xlabel(m1)
        ylabel(m2)
%         lims = quantile([ax.XLim ax.YLim], [0 1]);
%         xlim(lims)
%         ylim(lims)
        axis equal
        legend off

        F.print(F.prefix_better(sprintf('%s_%s', m1, m2)));
    end
    
    function [xx, yy, ax] = allP_directionality(F, metrics)
        % Shows the mean directionality index over 5s intervals
        
        data = F.CompiledData;
        if nargin < 2, metrics = F.Metrics; end
        
        figure('units', 'inches', 'position', [0 0 2 2]);
        mask = contains(data.metric, metrics);
        data = data(mask, :);

        m1 = metrics{1}; m2 = metrics{2};

        [G, pat, sz] = findgroups(data.patient, data.seizure);

        mask =@(mm) strcmpi(data.metric, mm);

        xx = splitapply(@nanmean, data.dir_ind_05(mask(m1)), G(mask(m1)));
        yy = splitapply(@nanmean, data.dir_ind_05(mask(m2)), G(mask(m2)));
        F.gscatter_pat(xx, yy, pat);
        
        disagree = F.Disagree;
        dd = cellfun(@(x) find(strcmpi(x, compose('%s %d', pat, sz))), disagree);
        
        ax = gca;
        hold(ax, 'on')
        plot(ax, xx(dd), yy(dd), 'k*', 'markersize', 6);
        hold(ax, 'off')
        
        title('Directionality index');
        xlabel(m1)
        ylabel(m2)
        axis square
        legend off

        F.print(F.prefix_better([m1 '_' m2]));
    end
    
    function [yy, ax] = allP_large_shift_times(F, metric, thresh)
        % [N, ax] = F.large_shift_times(metric='M10', thresh=pi/2)
        % Shows the number of large shifts for each patient
        if nargin < 2 || isempty(metric), metric = 'M10'; end
        if nargin < 3, thresh = pi/2; end
        metric = validatestring(metric, F.Metrics);
        
        data = F.CompiledData;
        mask = strcmpi(data.metric, metric) ...
            & data.dtheta_dt > thresh;
        data = data(mask, :);
        
        yy = data.time - data.iw_center;
        ax = F.gscatter_pat(yy, data.patient);
        title(sprintf('Large shifts (%s)', metric))
        ylabel('Time [s]')
        
        
    end
    
    function [xx, yy, ax] = allP_large_shift_counts(F, metrics, thresh)
        % [N, ax] = F.large_shift_counts(metric='M10', thresh=pi/2)
        % Shows the number of large shifts for each patient
        if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
        if nargin < 3, thresh = pi/2; end
        if numel(cellstr(metrics)) < 2; metrics = {metrics, ''}; end
        
        data = F.CompiledData;
        m1 = metrics{1}; m2 = metrics{2};
        mask = contains(data.metric, metrics);

        data = data(mask, :);
        
        [G, pat, ~] = findgroups(data.patient, data.seizure);
        
        mask =@(mm) strcmpi(data.metric, mm);
        
        if isempty(m2)
            xx = [];
            yy = splitapply(@(x) sum(x >= thresh), data.dtheta_dt(mask(m1)), G(mask(m1)));
            jitter = @(xx) xx + .3*(rand(size(xx)) - .5);
            ax = F.gscatter_pat(jitter(yy), pat);
            ylabel([m1 ' [counts]'])
        else
            xx = splitapply(@(x) sum(x >= thresh), data.dtheta_dt(mask(m1)), G(mask(m1)));
            yy = splitapply(@(x) sum(x >= thresh), data.dtheta_dt(mask(m2)), G(mask(m2)));

            jitter = @(xx) xx + .3*(rand(size(xx)) - .5);
            ax = F.gscatter_pat(jitter(xx), jitter(yy), pat);
            xlabel([m1 ' [counts]'])
            ylabel([m2 ' [counts]'])
        end
        
        title('Large shifts')
        axis square
        
        F.print(F.prefix_better([m1 '_' m2]));
        
    end
    
    function ln = allP_first_discharge_time(F, metrics)
        % Shows when the first discharge occurs relative to the IW (iw
        % center)
        data = F.CompiledData;
        T = F.get_iw_table;
        T = T(T.wave_num == T.main_wave, :);
        
        
        if nargin < 2, metrics = 'M'; end
        if numel(cellstr(metrics)) == 1, metrics = {'', metrics}; end
        
        mask = ismember(data.metric, metrics) & isfinite(data.iw_center);  % & isfinite(data.dir) 
        data = data(mask, :);
        
        [G, pat, sz, mm] = findgroups(data.patient, data.seizure, data.metric);
        
        t_val = splitapply(@min, (data.time - 1*data.iw_center), G);
        sp_val = nan(size(t_val));
        
        for ii = 1:length(pat)
            mask = strcmpi(T.patient, pat(ii)) & T.seizure == sz(ii);
            if T.nchannels(mask) < 30, continue; end
            sp_val(ii) = T.prop_sp(mask);
            
        end
        
        val = t_val .* sp_val;
        val = t_val;
        if numel(unique(mm)) == 1
            xx = [];
            yy = val;
        else
            pat = pat(strcmp(mm, metrics{1}));
            xx = val(strcmp(mm, metrics{1}));
            yy = val(strcmp(mm, metrics{2}));
        end

        ln = F.gscatter_pat(xx, yy, pat);
        legend('location', 'eastoutside')
        if ~isempty(metrics{1})
            xlabel(metrics{1})
            ylabel(metrics{2});
            title('First discharge time [s]')
%             drawnow;
            axis(gca, 'equal')
        else
            ylabel('[s]')
            title(sprintf('First discharge time (%s)', metrics{2}))
        end
        
    end
    
    
    function [yy, ax] = allP_iw_stats(F, min_electrodes)
        % Show info about detected IW for each patient and seizure
        
        if nargin < 2, min_electrodes = 0; end
        PCUTOFF = -8;
        
        fields = ["nchannels", "fr_peak_mu", "prop_sp", "p"];        
        h = figure('name', 'iw_stats', 'position', [0 0 .6*(numel(fields)+1) 1.5] * 2);

        data = F.get_iw_table(0);  % load the iw stats
        data.p = log10(data.p/5);  % rescale p-value so that it's in terms of p=5eXX
        
        % rescale very low pvalues
        lo_p = data.p < PCUTOFF;
        data.p(lo_p) = rescale(data.p(lo_p), PCUTOFF-1, PCUTOFF);

        % rescale low channel counts
        lo_ch = data.nchannels < min_electrodes;
        data.nchannels(lo_ch) = ...
            rescale(data.nchannels(lo_ch), min_electrodes-10, min_electrodes - 1);


        labels = struct(...
            'fr_peak_mu', {{{'Peak'; 'firing rate'}, '[Hz]', 'linear'}}, ...
            'nchannels', {{'Nchan', '[#]', 'linear', (min_electrodes:10:100)}}, ...
            'prop_sp', {{{'Speed'; '(log scale)'}, '[mm/s]', 'log', [0:.2:1, 2:2:4]}}, ...
            'moransI', {{{'Moran''s'; 'index'}, '[unitless]', 'linear'}}, ...
            'p', {{'p-val', 'p=5eYY', 'linear', PCUTOFF:0}}, ...
            'crossing_time', {{'Crossing time', '[s]', 'linear'}});

        data.disagree = ismember(compose('%s %d', string(data.patient), data.seizure), F.Disagree);

        T = tiledlayout(gcf, 1, numel(fields));
        for ff = fields
            ax = nexttile(T);
            if ff == "nchannels"
                dat = data; 
            else
                dat = data(data.nchannels >= min_electrodes, :); 
            end

            % show all waves in outline
            ln = F.gscatter_pat(dat.(ff), dat.patient);
            set(ln, 'markerfacecolor', 'none');
            hold(ax, 'on')
            for ll = ln', ll.ZData = -1*ones(size(ll.XData)); end


            % show disagrees with asterisks
            mask = dat.disagree & dat.wave_num == dat.main_wave;
            ln = F.gscatter_pat(dat.(ff)(mask), dat.patient(mask));
            set(ln, 'marker', '*', 'color', [0 0 0], 'markersize', 6);
            for ll = ln', ll.ZData = 1*ones(size(ll.XData)); end

            % show main wave in color
            [G, pat] = findgroups(dat.patient, dat.seizure);
            yy = splitapply(@(a, b, c) c(a == b), dat.wave_num, dat.main_wave, dat.(ff), G);
            ln = F.gscatter_pat(yy, pat);             
            for ll = ln', ll.ZData = 0*ones(size(ll.XData)); end
            
            if ff == "nchannels"
                yline(min_electrodes, 'r--');
            end
            if ff == "p"
                yline(PCUTOFF, 'k:');
            end

            if numel(labels.(ff)) == 4, yticks(labels.(ff){4}); end

            hold off
            set(ax, 'yscale', labels.(ff){3});
            set(findobj(T, 'type', 'line'), 'linewidth', 1)
            title(labels.(ff){1})
            ylabel(labels.(ff){2})
        end
        lgd = findobj(T, 'type', 'legend');
        set(lgd, 'visible', 'off');
        set(lgd(1), 'visible', 'on', 'location', 'eastoutside')
        lgd(1).String(contains(lgd(1).String, 'data')) = [];
        F.print(h, F.prefix_better(''));
    end
    
    function compare_metrics(F, pat, metrics)
        % Plot direction rasters from metrics for patient in pat
        if nargin < 2 || isempty(pat), pat = "c7_Seizure1"; end
        if isnumeric(pat), pat = string(F.SeizureInfo.name{pat}); end
        if nargin < 3, metrics = string(F.Metrics(1:2)); end
        
        sz_num = find(strcmpi(pat, F.SeizureInfo.name));
        
        h = figure('name', 'compare_fits', 'position', [0 0 2 1] * 2);
        T = tiledlayout(h, 3, 6, 'TileSpacing', 'compact');
        
        if pat == "c7_Seizure1", rotate_by = -pi/2; else, rotate_by = 0; end
            
        S1 = [1 3];
        S2 = [3 3];
        
        m1 = metrics(1);
        m2 = metrics(2);
        
        fit1 = WaveProp.load(F.get_fname(sz_num));
        fit2 = fit1.(m2);
        fit1 = fit1.(m1);
        
        fit1.RotateBy = rotate_by;
        fit2.RotateBy = rotate_by;
        
        % Density plot
        ax = nexttile(T, 4, S2);
        [f1, xi] = circ_ksdens(fit1.Direction);
        f2 = circ_ksdens(fit2.Direction, xi);
        
        plot(ax, rad2deg(xi), f1, 'displayname', m1, 'color', F.Style.(m1).color);
        hold(ax, 'on'); 
        plot(ax, rad2deg(xi), f2, 'displayname', m2, 'color', F.Style.(m2).color); 
        hold(ax, 'off');
        legend(ax);
        xticks(ax, -180:90:180);
        xlabel(ax, 'Direction [\circ]')
        ylabel(ax, 'PDF');
        grid(ax, 'on');
        
                
        % Direction v. time plots
        ax1 = nexttile(T, S1);
        fit1.direction_raster(ax1);
        ylabel(m1)
        xticklabels(ax1, []);
        title(ax1, '');
        
        ax2 = nexttile(T, S1);
        fit2.direction_raster(ax2);
        ylabel(m2);
        xticklabels(ax2, []);
        title(ax2, '')
        
        ax3 = nexttile(T, S1);
        if numel(fit1.time) < numel(fit2.time)
            fit1.direction_raster(ax3, fit1.diff(fit2));
        else
            fit2.direction_raster(ax3, fit2.diff(fit1));
        end
        ylabel('Diff');
        xlabel('Time [s]');
        title(ax3, '');
        
        aa = [ax1 ax2 ax3];
        linkaxes(aa, 'x');
        set(aa, 'ytick', [-180 0 180], 'yaxislocation', 'left')
        
        
        % Labels
        title(T, F.SeizureInfo.display_names{sz_num}, 'fontsize', 11);
        
        % Make a colorwheel
        pax = fit1.colorwheel;
        
        % Print
        pat = F.SeizureInfo.patient{sz_num};
        tag = sprintf('%s_%s_%s', pat, m1, m2);
        F.print(h, F.prefix_better(tag));
        F.print(pax.Parent, F.prefix_better('colorwheel'));
        close(pax.Parent);
        
    end
    
    function data = compile_fits(F, interval, dir_index_win)
        % int is the interval surrounding the IW to assign phase 1 (early)
        % or phase 3 (late)
        
        if nargin < 2 || isempty(interval), interval = Inf; end  
        if nargin < 3 || isempty(dir_index_win), dir_index_win = 5; end
        USE_MID = true;  % compute the intervals from the midpoint of the IW (or from the bounds)
        
        % Load the IW and TW data
        W = WaveProp.load([], F.Metrics);
        sz = F.SeizureInfo;
        data = table;
        iw_info_all = BVNY.get_iw_main;
        
        % In addition to metrics, compute "sub20"s (recompute TW using only
        % first 20 electrodes)
        % Was trying to show that TW propagate nonlinearly (by showing that
        % the initial propagation direction doesn't always match the
        % direction estimate using the full set of electrodes), but Mark thinks
        % this is obvious enough that this analysis is not necessary
%         metrics = [F.Metrics, ...
%             arrayfun(@(mm) sprintf("%ssub20", mm), F.Metrics)];
        metrics = F.Metrics;
        
        
        

        for mm = metrics  % loop through metrics
%             for ii = 1:height(sz)
            for name = string(sz.name)'
%                 name = sprintf('%s_Seizure%d', sz.patient{ii}, sz.seizure(ii));
                ii = find(contains(sz.name, name));
                iiW = cellfun(@any, strfind(W.Name, name));
                if ismember(name, fields(iw_info_all))
                    iw_info = iw_info_all.(name);
                    iw_info.nchan = numfinite(iw_info.template.template(iw_info.main_wave, :));
                else
                    iw_info = struct( ...
                        'center', nan, ...
                        'range', [nan nan], ...
                        'phi', nan, ...
                        'phi_pval', nan, ...
                        'nchan', 0);
                end
                

                % Get the TW info
                ss = strsplit(mm, 'sub');
                fit = W.(ss(1))(iiW);
                fit.MinFinite = BVNY.MinFinite(fit);
                
                % If this is a subset fit (i.e. "Msub20"), refit the TOA
                if numel(ss) == 2
                    fit = copy(fit); % these are handles; if you don't copy and you put the "sub" names first (which you probably won't), you will overwrite the originals and need to reload
                    fit.TOA = fit.toa_from_first_N_electrodes(str2double(ss(2)));
                    fit.refit_data;
                    fit.MinFinite = 5;  % Lower this requirement
                end
                
                
                [dir, sp] = fit.discharge_directions;
                data_temp = table( ...
                    dir, ...
                    fit.time, ...
                    sp, ...
                    [nan; fit.dtheta_dt], ...
                    fit.directionality_index(1), ...
                    fit.directionality_index(2), ...
                    fit.directionality_index(5), ...
                    fit.directionality_index(10), ...
                    'VariableNames', {'dir', 'time', 'speed', ...
                    'dtheta_dt', ...
                    'dir_ind_01', ...
                    'dir_ind_02', ...
                    'dir_ind_05', ...
                    'dir_ind_10', ...
                    });
                N = height(data_temp);
    %             name = strsplit(W.Name{ii}, {'_', 'Seizure'});
    %             data_temp.patient = repmat(string(name{1}), N, 1);
    %             data_temp.seizure = repmat(string(name{2}), N, 1);
                data_temp.patient = repmat(string(sz.patient{ii}), N, 1);
                data_temp.seizure = repmat(sz.seizure(ii), N, 1);
                data_temp.metric = repmat(mm, N, 1);
                data_temp.ubo = repmat(sz.UBO(ii), N, 1);
                
                phase = zeros(N, 1);
                iw0 = iw_info.range(1); % iw start
                iwF = iw_info.range(2); % iw end
                iwM = iw_info.center; % iw midpoint
                
                if USE_MID
                    phase(data_temp.time > iw0 - interval ...
                        & data_temp.time < iwM) = 1;  % pre-IW
                    phase(data_temp.time > iw0 ...
                        & data_temp.time < iwF) = 2;  % during IW
                    phase(data_temp.time > iwM ...
                        & data_temp.time < iwF + interval) = 3;  % post-IW
                    
                else
                    phase(data_temp.time > iw0 - interval ...
                        & data_temp.time < iw0) = 1;  % pre-IW
                    phase(data_temp.time > iw0 ...
                        & data_temp.time < iwF) = 2;  % during IW
                    phase(data_temp.time > iwF ...
                        & data_temp.time < iwF + interval) = 3;  % post-IW
                
                end
                
                data_temp.phase = phase;
                data_temp.iw_center = repmat(iwM, N, 1);
                data_temp.iw_angle = repmat(iw_info.phi, N, 1);
                data_temp.angle_p_val = repmat(iw_info.phi_pval, N, 1);
                data_temp.weight = 1/N * ones(N, 1);
                data_temp.nchannels = repmat(iw_info.nchan, N, 1);

                
                % Compare TW templates to IW templates
                if ~isnan(iw_info.center)
                    [rho, pval, tt] = fit.correlation(iw_info.template.template);
                    [~, locs] = min(abs(tt - fit.time'), [], 2);
                
                    temp = nan(N, 1);
                    temp(locs) = rho(:, iw_info.main_wave);  
                    data_temp.rho = temp;
                
                    temp = nan(N, 1);
                    temp(locs) = pval(:, iw_info.main_wave);
                    data_temp.rho_pval = temp;

                    temp = nan(N, 1);
                    temp(locs) = tt - iw_info.template.time(iw_info.main_wave);
                    data_temp.rho_time_offset = temp;
                else
                    [data_temp.rho, data_temp.rho_pval, data_temp.rho_time_offset] = ...
                        deal(nan(N, 1));
                end

                
                
                data = [data; data_temp];

            end
        end

    end

    function f = get_fname(F, ii)
        f = sprintf('%s_Seizure%d', ... 
            F.Seizures.patient{ii}, F.Seizures.seizure(ii));
    end
    
    function f = get_file(F, ii)
        % Wrapper for get_fname (which used to have a tail on it)
       f = F.get_fname(ii);
    end
    
    function outname = prefix_better(F, tag)
        caller = dbstack(1);
        assignin('base', 'caller', caller);
        assert(isstring(tag) || ischar(tag));
        outname = sprintf('%sfig6/%s_%s_%d', F.Prefix, ...
            strrep(caller(1).name, '.', '_'), ...
            tag, F.Smoothing); 
    end
    function outname = prefix(F, tag)
        assert(isstring(tag) || ischar(tag));
        outname = sprintf('%sfig6/fig6_%s_%d', F.Prefix, tag, F.Smoothing); 
    end
    function speed_box_summary(F)
        h = figure('units', 'inches', 'position', [0 0 2.35 3]);
        data_all = F.CompiledData;
%         [G, pat, sz, mtc, ubo] = findgroups(...
%             data_all.patient, ...
%             data_all.seizure, ...
%             data_all.metric, ...
%             data_all.ubo);    %#ok<ASGLU>
        [G, pat, sz, mtc, ubo, phase] = findgroups(...
            data_all.patient, ...
            data_all.seizure, ...
            data_all.metric, ...
            data_all.ubo, ...
            4 - data_all.phase);    %#ok<ASGLU>
        [~, mm] = findgroups(mtc);
        cmap = nan(numel(mm), 3);
        for ii = 1:numel(mm)
            cmap(ii, :) = F.Style.(mm(ii)).color;
        end
        cmap = [cmap; .5 * [1 1 1]];
        
        mask = ~(strcmpi(ubo, 'o') | phase == 2);
        mn = splitapply(@nanmean, data_all.speed, G);    
        boxplot(mn(mask), categorical([mtc(mask) ubo(mask) phase(mask)]), ...
            'datalim', [0 min(1800, quantile(mn(mask), 1-1.5/sum(mask)))], ...
            'extrememode', 'compress', ...
            'boxstyle', 'filled', 'medianstyle', 'target', ...
            'colors', cmap([1 4 1 4 2 4 2 4 3 4 3 4], :), 'factorgap', 10)
        
        ax = gca;
        ax.YGrid = 'on';
        ax.TickLength = [0 0];
        for m = mm'
            m1 = strcmpi(mtc, m) & strcmpi(ubo, 'u');
            m2 = strcmpi(mtc, m) & strcmpi(ubo, 'b');
            [~, p] = ttest2(mn(m1), mn(m2), 'vartype', 'unequal');
            mdl1 = fitdist(mn(m1), 'norm');
            mdl2 = fitdist(mn(m2), 'norm');
            display(p, sprintf('T-stat (unequal variance) %s', m));
            display(mdl1, sprintf('Unimodal norm fit: %s', m))
            display(mdl2, sprintf('Bimodal norm fit: %s', m))
        end
        F.print(h, [F.Prefix 'fig6/fig6_speed_box_summ'])
    end
    function dat = get.CompiledData(F)
        if isempty(F.CompiledData)
            F.CompiledData = F.compile_fits;
        end
        dat = F.CompiledData;
    end
    
    function make(F, patients, style)
        if nargin < 3, style = []; end
        if nargin < 2 || isempty(patients), patients = string(unique(F.Seizures.patient))'; end
        
        
        %%% F2: METHODS AGREE IN DATA %%%
        % Here we show that methods mostly agree in the data and
        % subsequently reduce the dataset to seizures where we see a
        % candidate IW that recruits at least 40 electrodes 
        
        % Some stats on agreement between methods
        F.allP_pdf_distance;
        F.allP_directionality;
        F.allP_proportion_finite;
        
        % Make examples for one agree patient and one disagree patient
        F.compare_metrics('c7_Seizure1');
        F.compare_metrics('c3_Seizure1');
        
        % only show stats for seizures where IW was detected on at least 30
        % electrodes. There are some cases where an IW is detected on
        % fewer electrodes but this is less likely to actually be an IW
        MIN_ELEC = 40;
        F.allP_iw_stats(MIN_ELEC); 
        

        
        %%% F3: TW DIRECTION WRT IW DIRECTION %%%
        % Here we show that there is no clear evidence to support the IW
        % hypothesis, but also not enough against it to reject
        
        % make the black and white direction v. time with IW plots for each
        % patient
        for pat = patients
            fprintf('Starting %s\n', pat{:})
            try F.hist_figs(pat, style); 
            catch ME, warning('make failed in %s', pat{:}); disp(ME); disp(ME.stack);
            end
        end
        F.save;
        F.dir_dist_by_phase();  % no longer used, but still in Illustrator so keep updating
        F.make_polarhist;  % show the pre-/post-IW polar histograms for each patient and seizure
        
        % Make the stick plots and polar histograms for all patients
        % inf is the interval surrounding the IW; I chose 40
        % (min_electrodes) since the main wave of most seizures appears on
        % at least 50 (40/50 makes no difference)
        F.pre_post_direction_stick_plots(inf, MIN_ELEC);  
        F.pre_post_direction_hist(inf, MIN_ELEC);
        
        
        %%% F4: TW/IW CORRELATION %%%
        % Here we show that [...TBD]
        
        % Show correlation v. time (relative to IW)
        F.allP_corr_v_time('M', 5e-2, MIN_ELEC);
        
        % Test for times of significance (this one takes a while)
        F.allP_hist2d_tw_v_iw_effectsize('M', 5e-2, MIN_ELEC, 'rho');
        
        % Show how many TW have a strong relationship to the main IW. Only
        % show these numbers for IW that appear on more than 40 electrodes.
        F.allP_quantiles_pval_tw_v_iw('M', 5e-2, MIN_ELEC);
        
        % Show all the waves from c7
        for pp = ["c7_Seizure1", "MG49_Seizure43", "MG63_Seizure4", "CUCX3_Seizure6"]
            for cval = ["pos" "neg" "zero"]
                F.show_all_sig_corr(pp, 'M', cval);  % c7, why you look different? Because non-linearity
            end
        end
        
        % Compare the directions determined using the full time window
        % versus just the first 20 observed times (c7, for example, follows
        % a curved trajectory and the corr/anti-corr looks like it
        % disagrees with the direction results.
        
        % Was trying to show that TW propagate nonlinearly by showing that
        % the direction estimate using only the first 20 electrodes doesn't
        % always match the estimate using the full set of electrodes, but
        % Mark thinks this is obvious enough that we don't need to show it.
        % If you put this back in, you will need to uncomment the line
        % where you add the "sub20" metrics.
        if 0  
        for pp = patients %#ok<UNRCH>
            h = F.show_dir_v_time_early(pp, [], 'raster');
            close(h);
            % Not sure yet if this one or the previous is more informative
            h = F.show_dir_v_time_early(pp, [], 'scatter');
            close(h);
        end
        end
        
        
        
 

    end
    
    
    function save(F, patients)
        if nargin < 2, patients = fieldnames(F.H); end
        if ischar(patients), patients = {patients}; end
        for pat = patients', save_fig_(F, F.H.(pat{:})); end
    end
    function close(F, metrics)
        if nargin < 2, metrics = fieldnames(F.H)'; end
        for mm = string(metrics), close(F.H.(mm)); end
        F.H = [];
    end
    
    function [h, ax] = create_fig(F, patient)
        C = numel(F.Metrics) + 1;
        R = sum(strcmpi(F.Seizures.patient, patient));
        pos = F.subplot_position(R, C);
        h = figure( ...
            'units', 'inches', ...
            'position', [0, 0, ...
                sum(pos([1 3])) + 2*F.Margin(1), ...
                sum(pos([2 4])) + 2*F.Margin(2) ...
                ], ...
            'resize', 'off');
        
        ax = gobjects(R, C);
        for r = 1:R
            for c = 1:C
                ax(r, c) = axes(h, 'units', 'inches', 'ytick', [], ...
                    'position', [F.Margin 0 0] + F.subplot_position(r, c));
%                 title(ax(r, c), 'title');
            end
        end
        ax = flipud(ax);
        ylabel(axes(h, ...
                'visible', 'off', ...
                'units', 'norm', ...
                'position', [0 0 1 1], 'tag', 'label'), ...
            'Patient', ...
            'position', [0 .5 0], ...
            'visible', 'on', ...
            'verticalalignment', 'top', ...
            'tag', 'patient');
        
        
    end
    
    
    function pos = subplot_position(F, r, c)
        
        c = c-1; r = r-1; % zero-indexing
        
        w = 1.25;
        h = .5;
        
        pos = [ ...
            (w+2*F.Padding) * c + F.Padding, ...
            (h+2*F.Padding) * r + F.Padding, ...
            w, h];
        if c == numel(F.Metrics), pos(3) = .3; end
    end
    
    hist_figs2(F, patient)
    hist_figs(F, patient, style)
    hist_figs_speed(F, patient)
    dir_dist_summary(F, metrics, period)
    dir_dist_by_phase(F, metrics, class)
    plot_ttests(F)
end

methods (Static)
   [files, names] = txt2files(fname) 
   style = set_style
   
   function T = get_iw_table(min_electrodes)
       if nargin < 1, min_electrodes = 0; end
        T = readtable('iw_table');
        T(T.nchannels < min_electrodes, :) = [];
        
        T.seizure_dur = T.onset ./ T.onset_rel;
        T.phi_alt_deg = T.phi_alt/pi * 180;
        T.crossing_time = T.prop_sp;
        T.prop_sp = 10 * sqrt(2) * .4 ./ T.prop_sp;
        bad_speed = isoutlier(T.prop_sp, 'ThresholdFactor', 100);
        disp(find(bad_speed))
        T(bad_speed, :) = [];
        T.vmr_increase = T.vmr_median_iw ./ T.vmr_median_late;
        
        G = findgroups(T.patient, T.seizure);
        T.w_num = zeros(height(T), 1);
        for ii = 1:max(G)
            [~, T.w_num(G==ii)] = sort(T.onset(G==ii));
        end

   end
    
   function iw_info = iw_info_all()
        iw_info = foldercontents2struct('iw_info');
   end
    
   function [sc] = gscatter_pat(xx, yy, pat, ax)
        
        % only 6 because MG49 only shows up once now and this way colors
        % match with fig0
        cmap = repmat(lines(6), 2, 1);
        lso = 'oooooo^^^^^';
        
        switch nargin
            case 4
            case 3 
                if isa(pat, 'matlab.graphics.axis.Axes') %(F, yy, pat0, ax)
                    ax = pat;
                    pat = yy;
                    yy = xx;
                    xx = [];
                else  % (F, xx, yy, pat0);
                    ax = gca;
                end
            case 2 % (F, yy, pat0)
                pat = yy;
                yy = xx;
                xx = [];
                ax = gca;
            otherwise
                error('Expected F.(xx*, yy, pat, ax*). Asterisks denote optional.')
        end
        assert(iscell(pat) || isstring(pat));
        
        
        % Create mapping and reorder data
        sz = SeizureInfo;
        [~, p, pa] = findgroups(sz.patient, sz.patientAlt);
        MAP = containers.Map(p, pa);
        G2 = cellfun(@(x) MAP(x), pat);
        
        [G2, so] = sort(G2);
        yy = yy(so);
        pat = pat(so);
        
        % maintain the same colors/lines every time
        cmap = cmap(unique(G2), :);
        lso = lso(unique(G2));

        if isempty(xx)
            xx = rescale(G2, .8, 1.2, ...
                'inputmin', 1, 'inputmax', 11);
            XLIM = [.6 1.4];
            XTICK = [];
        else
            xx = xx(so);
            XLIM = [];
            XTICK = 'auto';
        end

        for ii = 1

            sc = gscatter(ax, xx, yy, pat, ...
                cmap, lso, [], 'on');
            for ss = sc', ss.MarkerFaceColor = ss.Color; end
            
            grid on
            box off
            xticks(XTICK);
            xlabel([]);
            if ~isempty(XLIM), xlim(XLIM); end
            
        end

    end

end
    
    
end


%% --- Local functions ---

function save_fig_(F, h)

outname = F.prefix(h.Tag);
F.print(h, outname); 
% savefig(h, outname)


end




