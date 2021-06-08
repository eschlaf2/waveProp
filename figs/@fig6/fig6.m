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
    SlidingDirs
    Disagree
    WaveFits
    ExampleSeizures = ["c7_Seizure1" "MG49_Seizure43" "BW09_Seizure2" "CUCX4_Seizure1"]
end

properties
    Margin = [0.25 0.25];
    Padding = 0.1;
end

methods
    function F = fig6
        BVNY.set_style();
        F.Seizures = F.SeizureInfo;
%         mask = contains(F.SeizureInfo.class, 'TLE');
%         F.Seizures = F.SeizureInfo(mask, :);
    end
    
    function T = show_corr_dir_mismatches(F, metric, min_electrodes)
        % Show TW patterns for outliers in the correlation-direction
        % relationship
        % T = show_corr_dir_mismatches(F, metric="M", min_electrodes=50)
        
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(min_electrodes), min_electrodes = 50; end
        metric = validatestring(metric, F.Metrics);
        
        
        % Parameters for figures
        MAKE_ONE_PAGE = true;  % only show the first page of discharges
        rows = 4; cols = 8;
        panes_per_fig = rows * cols;
        
        % Get data
        data = F.CompiledData;
        mask = isfinite(data.dir) ...  % Reduce to data with significant direction estimates
            & data.nchannels >= min_electrodes ...  % ... with strong IW
            & data.metric == metric;  % ... and matches <metric>
        data = data(mask, :);
        
        % Rename some fields; convert from radians to degrees
        data.dir0_d = rad2deg(data.dir);
        data.iw_angle_d = rad2deg(data.iw_angle);
        data.dir_d = rad2deg(fix_angle(data.dir-data.iw_angle));
        data.abs_dir_d = abs(data.dir_d);
        
        % Fit a regression and identify points with large residuals; 
        [~, ~, r] = regress(data.rho, ...  % get the residuals
            [ones(size(data.abs_dir_d)) data.abs_dir_d]);
        big_r = abs(normalize(r)) > 2;  % ... and identify outliers
        
        % Reduce data to the set with large, counterindicative residuals
        mask_ = big_r;  % next lines commented out bc of c7
%         & ( ...  % outlier from regression and 
%                 (abs(data.dir_d) > 90 & data.rho > 0) ...  % ... dir near 180 and corr
%                 | (abs(data.dir_d) < 90 & data.rho < 0) ...  % ... or dir near 0 and anti-corr
%             );
        
        
        data = data(mask_, :);
        data.resid = r(mask_);
        
        % Sort data by abs(resid)
        [~, so] = sort(abs(data.resid), 'descend');
        data = data(so, :);
        
        % If there are too many outliers, only show the first page
        num_found = height(data);
        
        if num_found > panes_per_fig && MAKE_ONE_PAGE
            G = findgroups(data.patient);
            
            Gcount = zeros(size(G));
            for ii = 1:length(Gcount)
                Gcount(ii) = sum(G(1:ii) == G(ii));
            end
            
            p_max = cols;  % no more than one row of examples from each patient.
            mask_ = Gcount <= p_max;
            while sum(mask_) < panes_per_fig
                p_max = p_max + 1;
                mask_ = Gcount <= p_max;
            end
            dat_sub = data(mask_, :);
            data = dat_sub(1:panes_per_fig, :);
            
        end
        
        
        % Get the WaveProp objects
        [G, pat, sz] = findgroups(data.patient, data.seizure);
        M = cell(numel(pat), 1);
        for ii = 1:length(M)
            fit_temp = F.get_wave_fit(sprintf('%s_Seizure%d', pat(ii), sz(ii)));
            M{ii} = fit_temp.(metric);
        end
        
        
        % Make the figures
        num_figs = ceil(height(data) / panes_per_fig);  % max (ros*cols) panes per figure
        for ii = 1:num_figs
            h(ii) = figure('name', sprintf('rho_dir_mismatches_%d', ii)); fullwidth(true); %#ok<AGROW>
            T(ii) = tiledlayout(h(ii), rows, cols); %#ok<AGROW> % If you use 'flow' here, the scatters will have be the wrong size
            title(T(ii), sprintf("Mismatched \\rho-\\phi (%s, %d/%d)", ...
                metric, height(data), num_found));
        end
        
        % Plot the panels
        for ii = 1:height(data)
            which_fig = ceil(ii / panes_per_fig);
            ax = nexttile(T(which_fig), ii); 
            M{G(ii)}.plot2D(ax, data.time(ii)); 
            hold(ax, 'on');  % add an invisible dot to make a pretty annotation
            ln = plot(ax, 0, 0, '.', 'color', 'none'); 
            hold(ax, 'off');
            title(ax, sprintf("%s %d", pat(G(ii)), sz(G(ii))));
            ant_str = sprintf( ...
                '\\rho:  %3.2f,  \\phi:  %3.0f°\n\\phi TW (IW):  %3.0f° (%3.0f°)\n', ...
                data.rho(ii), data.dir_d(ii), data.dir0_d(ii), data.iw_angle_d(ii));
            legend(ln, ant_str, 'location', 'southoutside', ...
                'color', 'none', 'box', 'off');
            
            box(ax, 'on');
            axis(ax, 'square')
            xlim(ax, [0 11])
            ylim(ax, [0 11])
            
        end
        
        % Match color scheme to method
        all_ax = findobj(T, 'type', 'axes');
        len = size(all_ax(1).Colormap, 1);
        mtc_color = F.Style.(metric).color;
        cmap = make_diverging_colormap(...
            [brighten(mtc_color, .9); brighten(mtc_color, -.9)], ...
            mtc_color, len);
        set(all_ax, 'colormap', cmap);
        
        
        % If there aren't enough examples to fill the rows of T, copy an axis to the
        % end. Otherwise, the figure will change size and screw up your
        % clips in Illustrator.
        if ii < (rows - 1) * cols + 1
            ax = nexttile(T, panes_per_fig);
            M{G(ii)}.plot2D(ax, data.time(ii)); 
            % I don't think any of the prettify above affects the sizing so
            % nothing else should be needed.
            set(ax, 'colormap', gray);  % setting to gray just for visibility
        end
        
        
        % Match the size of the scatters for each plot
        sc = findobj(T, 'type', 'scatter');
        set(sc, 'sizedata', M{1}.get_dot_size(ax));
        
        % Fix legend locations
        % You can't move legends inside tiledlayouts. I feel like this is
        % usually more straightforward, but right now this is what it is...
        
        % Copy children of T into figure
        kiddos = T.Children;  
        pos0 = cell2mat(arrayfun(@(ax) ax.Position, kiddos, 'uni', 0));
        new_kids = copyobj(T.Children, h);
        for ii = 1:numel(new_kids)
            new_kids(ii).Position = pos0(ii, :);
        end
        delete(kiddos);
        
        % Update legend locations and make font size a little bigger
        align_right = @(ax, lgd) ...  % Get the "left" legend position
            sum(ax.Position([1 3])) - lgd.Position(3);
        align_top = @(ax, lgd) ...  % Get the "bottom" legend position
            ax.Position(2) - lgd.Position(4);
        for aa = findobj(h, 'type', 'axes')'
            lgd = aa.Legend;
            lgd.Position(1) = align_right(aa, lgd);
            lgd.Position(2) = align_top(aa, lgd);
            lgd.FontSize = 7.5;
        end
        
        % dunno why, but if you don't `drawnow` the following warning pops
        % up (although export actually looked ok)
        %   Warning: Exported image displays axes toolbar. 
        %   To remove axes toolbar from image, export again. 
        drawnow;  
        
        for ii = 1:num_figs
            tag = sprintf('%s_%d', metric, ii);
            F.print(h(ii), F.prefix_better(tag));
        end
    end

    function T = show_all_sig_corr(F, sz_name, metric, pos_or_neg, time_range)
        % Show each significant wavefit/correlation pair for a given
        % seizure number and metric
        if nargin < 2 || isempty(sz_name), sz_name = 'c7_Seizure1'; end
        if nargin < 3 || isempty(metric), metric = 'M'; end
        if nargin < 4 || isempty(pos_or_neg), pos_or_neg = "pos"; end
        if nargin < 5 || isempty(time_range), time_range = [-inf inf]; end
        
        
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
        M = F.get_wave_fit(F.get_file(sz_num)); M = M.(metric);
        [iw, mw] = BVNY.get_iw_info(sz_num);
        tpl = iw{mw}.template.template(mw, :, :);
        
        
        % get the correlations and directions (only show sig fits for both)
        [rho, rho_p, rho_t] = M.correlation(tpl);
        
        % Get only those correlations in the specified time range
        tr = time_range + iw{mw}.center;
        mask = rho_t >= tr(1) & rho_t <= tr(2); 
        rho = rho(mask); rho_p = rho_p(mask); rho_t = rho_t(mask);
        
        
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
            fit = F.get_wave_fit(fname);
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

                out.(mm).pre(gg, :) = [circ_mean(pre) numel(pre)];
                out.(mm).post(gg, :) = [circ_mean(post) numel(post)];
                out.iw_angle_pval(gg) = dat.angle_p_val(1);
            end
        end
        
    end
    
    
    function out = kl_divergence(F)
        % This computes a bunch of measures comparing the methods and the
        % resulting distributions; not just kl_divergence. 
        % Used to generate scatter plots
        
        % Dkstat: This shows a circular variation of the KS statistic
        % comparing the distribution of the differences to a Von Mises
        % distribution with mean and kappa matching the sample mean and
        % kappa (i.e. circular equivalent to goodness-of-fit between
        % differences and normal distribution). Lower Kstat means the
        % differences are better described by the VM.
        
        RES = 64;
        XI = linspace(-pi, pi, RES + 1);
        ANGULAR_LAGS = diff(XI(1:2)) .* (-numel(XI)+1:numel(XI));
        MTC = F.Metrics;
        WIN = 10;  % number of seconds in sliding window
        QQ = [.05 .95];  % quantiles to use for moving stat summaries
        
        data = F.CompiledData;
        mask = ismember(data.metric, MTC);
        data = data(mask, :);
        
        [G, pat, sz] = findgroups(data.patient, data.seizure);
        

        % Functions to compute Kullback-Leibler divergence and Kuiper
        % statistic
        kl_dist = @(p, q) nansum(p .* log(p./q));  % units in nats
        kstat = @(p, q) (max([cumsum(p(:)) - cumsum(q(:)); 0]) ...
            + max([cumsum(q(:)) - cumsum(p(:)); 0]));
        
        W = F.WaveFits;
        mpairs = nchoosek(MTC, 2);
        nP = size(mpairs, 1);
        
        [Dkl, Dmode, Dmean, ...
            mean_diffs, std_diffs, skew_diffs, kappa_diffs, ...
            var_diffs, mov_var_diffs, ...
            mov_mean_diffs, mov_kstat, ...
            corr_pval, corr_coeff, kstat_diffs, pdf_kstat, mutual_info, ...
            pdf_xcorr_coeff_normed, pdf_kstat_corr_aligned, ...
            pdf_xcorr_coeff, pdf_xcorr_angular_lag] = deal(nan(max(G), nP));
        
        [bias, variance, CC, mov_std_diffs] = deal(nan(max(G), 3, nP));
        
        for ii = 1:max(G)  % for each seizure
            
            % Compute KL divergence of diffs from VM normal
            % Compare dist of diffs to VM normal with same center and std
%             W = WaveProp.load(sprintf('%s_Seizure%d', pat{ii}, sz(ii)), F.Metrics);
            sz_num = find(strcmpi(W.Name, sprintf('%s_Seizure%d', pat(ii), sz(ii))));
            
            
            for jj = 1:nP  % For each pair of metrics
                
                
                % Get the right fit objects
                fit1 = W.(mpairs{jj, 1})(sz_num);
                fit2 = W.(mpairs{jj, 2})(sz_num);
                                            
                %%% Compute stats about the distributions
                pdf1 = circ_ksdens(fit1.Direction, XI);
                pdf2 = circ_ksdens(fit2.Direction, XI);
                
                % normalize so they sum to 1
                pdf1 = pdf1/sum(pdf1); 
                pdf2 = pdf2/sum(pdf2);
            
                [~, loc] = max([pdf1; pdf2], [], 2);
                modes = XI(loc);
                Dmode(ii, jj) = angle(exp(1j* diff(modes) ));  % return the angle between -pi and pi
                means = cellfun(@(s) circ_mean(rmmissing(s.Direction)), {fit1 fit2});
                Dmean(ii, jj) = angle(exp(1j* diff(means) ));
            
                % compute the KS stat between the distributions
                pdf_kstat(ii, jj) = kstat(pdf1, pdf2);
                
                % Compute the circular cross correlation and determine the
                % shift
                cc = cconv(pdf1, fliplr(pdf2), numel(ANGULAR_LAGS));
                [cc_max, max_ind] = max(cc);
                angular_lag = ANGULAR_LAGS(max_ind);
                pdf_xcorr_coeff(ii, jj) = cc_max;
                pdf_xcorr_angular_lag(ii, jj) = angle(exp(1j*angular_lag));  % return the angle between -pi and pi
                pdf_xcorr_coeff_normed(ii, jj) = max(normalize(cc));
                shift_by = round(angular_lag / diff(XI(1:2)));
                pdf_kstat_corr_aligned(ii, jj) = kstat(pdf1, circshift(pdf2, shift_by));
                
                %%% Compute stats about the per-discharge differences
                if numel(fit1.time) < numel(fit2.time) % note that this will resample the times to match the method with fewer sample points
                    diffs = fit1.diff(fit2);
                    diff_times = fit1.time;
                    
                    % Get time resampled directions
                    dir1 = fit1.Direction;
                    fit2.resample_t0(fit1.time);
                    dir2 = fit2.Direction;
                    fit2.Inds = [];  % set back as they were
                else
                    diffs = fit2.diff(fit1);
                    diff_times = fit2.time;
                    
                    % Get time resampled directions
                    dir2 = fit2.Direction;
                    fit1.resample_t0(fit2.time);
                    dir1 = fit1.Direction;
                    fit1.Inds = [];  % set back as they were
                end
                dir1(isnan(diffs)) = [];
                dir2(isnan(diffs)) = [];
                diff_times(isnan(diffs)) = [];
                diffs(isnan(diffs)) = [];
                
                % Compute mutual information
                edges = linspace(-pi, pi, 65);
                px = histcounts(fit1.Direction, edges, 'Normalization', 'prob');
                py = histcounts(fit2.Direction, edges, 'normalization', 'prob');
                pxy = histcounts2(dir1, dir2, edges, edges, 'normalization', 'prob');
                mutual_info(ii, jj) = sum(kl_dist(pxy, px'*py));
                
                
                % Compute stats on the differences between the fits
                % (sliding window)
                test_times = diff_times(1):1:diff_times(end);  % compute stats at 1 s intervals throughout the seizure
                
                t_win = @(tt) abs(diff_times - tt) <= WIN/2;  % Get a 10 s window of diffs
                N = arrayfun(@(tt) sum(t_win(tt)), test_times);
                lowN_mask = N < WIN;  % at least one discharge/sec
                test_times(lowN_mask) = [];
                mov_stats = arrayfun(@(tt) circ_stats(diffs(t_win(tt))), test_times, 'uni', 0);
                mov_cc = arrayfun(@(tt) circ_corrcc(dir1(t_win(tt)), dir2(t_win(tt))), test_times);
                
                CC(ii, :, jj) = quantile(mov_cc, [.5 QQ]);
                
                mov_mn = cellfun(@(s) s.mean, mov_stats);  % Extract the mean from each interval centered at test_time
                [~, pk] = max(abs(mov_mn));  % Get the peak bias in the seizure (not used currently; <bias> used instead)
                mov_mean_diffs(ii, jj) = mov_mn(pk);
                bias(ii, :, jj) = quantile(abs(mov_mn), [.5 QQ]);
                
                mov_std = cellfun(@(s) s.std, mov_stats);
%                 [~, pk] = max(mov_std);
                mov_std_diffs(ii, :, jj) = quantile(mov_std, [.5 QQ]);
                
                mov_var = cellfun(@(s) s.var, mov_stats);  % Extract the variance at each test_time; You may be tempted to use 2*var to convert from circular to angular variance, but don't. If you convert to degrees, double the variance then, but just looking at the variance is more interpretable on a scale of 0 to 1.
                [~, pk] = max(mov_var);
                mov_var_diffs(ii, jj) = mov_var(pk);
                variance(ii, :, jj) = quantile(mov_var, [.5 QQ]);
                
                
                
                circ_pdf = @(tt) circ_ksdens(diffs(t_win(tt)), XI);
                sample_pdfs = arrayfun(@(tt) circ_pdf(tt), test_times, 'uni', 0);
                sample_pdfs = cellfun(@(s) s / sum(s), sample_pdfs, 'uni', 0);  % normalize
                mov_kappa = arrayfun(@(tt) circ_kappa(diffs(t_win(tt))), test_times);
                vm_pdfs = arrayfun(@(mn, kappa) ...
                    circ_vmpdf(XI, mn, kappa), mov_mn, mov_kappa, 'uni', 0);
                vm_pdfs = cellfun(@(s) s/sum(s), vm_pdfs, 'uni', 0);
                mov_k = cellfun(@(p, q) kstat(p, q), sample_pdfs, vm_pdfs);
                mov_kstat(ii, jj) = nanmean(mov_k);
                
                
                % Compute stats for full distribution
                stats = circ_stats(diffs);
                kappa_diffs(ii, jj) = circ_kappa(diffs);
                vm_pdf = circ_vmpdf(XI, stats.mean, kappa_diffs(ii, jj));  % kappa);
%                 vm_pdf = circ_vmpdf(xi, 0, 1);  % kappa);
                vm_pdf = vm_pdf / sum(vm_pdf);  % normalize to sum to 1
                pdfD = circ_ksdens(diffs, XI)';
                pdfD = pdfD / sum(pdfD);  % normalize to sum to 1
                
                Dkl(ii, jj) = kl_dist(pdfD, vm_pdf);
                mean_diffs(ii, jj) = stats.mean;
                var_diffs(ii, jj) = stats.var;
                std_diffs(ii, jj) = stats.std;
                skew_diffs(ii, jj) = stats.skewness;
                kstat_diffs(ii, jj) = kstat(pdfD, vm_pdf);  % ...  % Kuiper stat
                
                % Compute circular corr coeff
                [cc, pval] = circ_corrcc(dir1, dir2);  % this is computing the 
                corr_coeff(ii, jj) = cc;
                corr_pval(ii, jj) = pval;
            end
            
            
            
        end
        
        out = table(Dkl, Dmode, Dmean, mean_diffs, std_diffs, var_diffs, skew_diffs, ...
            bias, variance, CC, ...
            mov_mean_diffs, mov_std_diffs, mov_var_diffs, mov_kstat, ...
            corr_pval, corr_coeff, kappa_diffs, kstat_diffs, pat, sz, ...
            mutual_info, pdf_kstat, pdf_xcorr_coeff, pdf_xcorr_angular_lag, ...
            pdf_kstat_corr_aligned, pdf_xcorr_coeff_normed);
        
    end
    
    function W = get.WaveFits(F)
        if isempty(F.WaveFits)
            W = WaveProp.load([], F.Metrics);
            F.WaveFits = W;
        end
        W = F.WaveFits;
    end
    
    function out = get.FitDiffs(F)
        if isempty(F.FitDiffs)
            W = warning; warning off;
            out = F.kl_divergence;
            out.Dmode_d = out.Dmode / pi * 180;
            out.pdf_shift_d = rad2deg(out.pdf_xcorr_angular_lag);
            out.mean_diffs_d = rad2deg(out.mean_diffs);
            
%             [mn, ul, ll] = circ_mean(out.mean_diffs);
%             MMM = rad2deg([mn, ll, ul]);
            
            
            
            F.FitDiffs = out;
            warning(W);
        end
        out = F.FitDiffs;
    end
    function W = get_wave_fit(F, name)
        ind = strcmpi(F.WaveFits.Name, string(name));
        for ff = fieldnames(F.WaveFits)'
            W.(ff{:}) = F.WaveFits.(ff{:})(ind);
        end
        
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
            W = F.get_wave_fit(F.SeizureInfo.name(ii));

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
        
        % For each phase show a polar histogram
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
        
            pax.ThetaAxisUnits = 'deg';
            lbl = pax.ThetaTickLabel;
            lbl(~ismember(lbl, {'0', '180', '90', '270'})) = {''};
            lbl = strrep(lbl, '270', '-90');
            lbl = strrep(lbl, '0', '0\circ');
            pax.ThetaTickLabel = lbl;
            pax.LineWidth = 1;        
            title(phase)
        end

        % Summarize results
        num_seiz = numel(out.patient);
        % definitions
        chg_fun = @(diff) (abs(diff) > pi/2);
        align_fun = @(pre) (abs(pre) < pi/2);
        iw_fun = @(pre, diff) (chg_fun(diff) & align_fun(pre));
        % apply definitions to both TW methods
        supportIW = arrayfun(@(mm) sum(iw_fun(out.(mm).pre(:, 1), out.(mm).diff)), F.Metrics);
        big_change = arrayfun(@(mm) sum(chg_fun(out.(mm).diff(:, 1))), F.Metrics);
        iw_aligned = arrayfun(@(mm) sum(align_fun(out.(mm).pre(:, 1))), F.Metrics);
        fprintf('Early TW-IW alignment (pre-IW < 90°): %d[%d]/%d seizures, %s[%s]\n', ...
            iw_aligned, num_seiz, F.Metrics);
        fprintf('Change in direction > 90°: %d[%d]/%d seizures, %s[%s]\n', ...
            big_change, num_seiz, F.Metrics);
        fprintf('Support IW-hyp (pre-IW < 90° & change > 90°): %d[%d]/%d seizures, %s[%s]\n', ...
            supportIW, num_seiz, F.Metrics);
        
        
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

        cols = @(mm) F.Style.(mm).color;
        mask = out.iw_angle_pval < .05;

        sz = SeizureInfo;
        [~, p, pa] = findgroups(sz.patient, sz.patientAlt);
        MAP = containers.Map(p, pa);
        G = cellfun(@(x) MAP(x), out.patient);
        [~, so] = unique(G);
        
        
        Gx = findgroups(G);  % use findgroups to remove empty patient groups
        for mm = string(F.Metrics)

            % make x-values out of patient number and seizure index
            xx = [Gx Gx]' + .25*[-1; 1];  


            for offset = [-360 0 360]
                % unwrap differences and convert to degrees
                yy = unwrap([out.(mm).pre(:, 1) out.(mm).post(:, 1)]') / pi * 180;  
                yy = yy + offset;
                N = [out.(mm).pre(:, 2) out.(mm).post(:, 2)]';
                yy(N == 0) = nan;
%                 yy = [yy [0; 180]]; %#ok<AGROW>


                % Indicate IW and anti-IW direction
                ls = ':-:';
                theta = [-180 0 180];
                for ii = 1:numel(theta)
                    yline(theta(ii), ['k' ls(ii)], 'linewidth', 1);
                end
                hold on
                plot(ax, xx, yy, 'color', cols(mm), ...
                    'marker', '.', 'markersize', 15, ...
                    'displayname', mm)


                % Indicate good IW direction fits
                gray_ = 0.1 * [1 1 1];
                plot(xx(1, mask), yy(1, mask), '.', 'color', gray_);
                plot(xx(2, mask), yy(2, mask), '.', 'color', gray_);

            end
        
        end
        
       
        % Prettify
        xticks(ax, 1:max(Gx) + 1)
        xticklabels(ax, [out.patient(so); 'HYP'])
        xtickangle(ax, 45)
        yticks(ax, -360:90:360)
        yticklabels(ax, num2str([0 90 -180 -90 0 90 180 -90 0]'))
        ylim([-1 1] * 200);
        ylabel('Mean direction (\circ)')
        grid on

        ln = findobj(ax, 'type', 'line');
        ln1 = findobj(ln, 'displayname', F.Metrics(1));
        ln2 = findobj(ln, 'displayname', F.Metrics(2));
        
         % Add the hypothesis
        xHyp = max(xx, [], 2) + 1;
        lHyp1 = copyobj(ln1(1), ax);
        set(lHyp1, 'xdata', xHyp, 'ydata', [0 180], 'color', .35*[1 1 1], ...
            'displayname', 'IW Hyp');
        lHyp2 = copyobj(ln1(1), ax);
        set(lHyp2, 'xdata', xHyp, 'ydata', [0 -180], ...
            'linestyle', '--', 'color', .5*[1 1 1]);
        
        % Add the legend
        legend([ln1(1) ln2(1) lHyp1], ... % string(F.Metrics), ...
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
    
    function sub = FitDiffs_sub(F)
        out = F.FitDiffs;
        sub = out(:, ["pat" "sz"]);
        sub.Dmean_d = rad2deg(out.Dmean);
        sub.pdf_kstat = out.pdf_kstat_corr_aligned;
        sub.bias_d = rad2deg(out.bias);
        sub.variance = out.variance;
        sub.xcorr_lag = rad2deg(out.pdf_xcorr_angular_lag);
        sub.xcorr_cc = out.pdf_xcorr_coeff;

    end
    
    function disagree = get.Disagree(F)
        out = F.FitDiffs_sub;
        
%         agree_fun = @() ( ...
%             abs(out.Dmean_d) < 60 ...
%             & out.pdf_kstat < .4 ...  % distributions match
%             );
        agree_fun = @() ( ...
            abs(out.xcorr_lag) < 60 ...
            & out.pdf_kstat < .4 ...  % distributions match
            );
%         agree_fun = @() (...
%             abs(out.xcorr_lag) < 60 ...
%             & out.xcorr_cc > 0 ...  % this is doing nothing...
%             );
        
        out.agree = agree_fun();
        out.Row = compose('%s %d', out.pat, out.sz);
        
        show_mads = ["Dmean_d" "pdf_kstat"];
        c=-1/(sqrt(2)*erfcinv(3/2));  % scale factor for scaled mad
        scaled_mad = @(dat) c * nanmedian(abs(dat - nanmedian(dat)));
        
        fprintf('Disagree:\n');
        fprintf('3 scaled MAD:\n');
        arrayfun(@(ff) fprintf('%15s: %0.3g\n', ff, 3*scaled_mad(out.(ff))), show_mads);
        disp(agree_fun)
        disagree = out(~out.agree, :);

    end
    
    function ax = show_distribution_comparison(F, ax, display_option, which_stat)
        % ax = show_distribution_comparison(F, ax=gca, display_option="corr", which_stat="ci")
        % Plots mean vs. kstat for each seizure
        % which_stat is fed to summarize_stat (use 'ci' to show CI, 'std'
        % to show STD, 'range' to show median and range)
        
        if nargin < 2, ax = gca; end
        if nargin < 3 || isempty(display_option), display_option = "corr"; end
        if nargin < 4, which_stat = []; end
        if ~isa(ax, 'matlab.graphics.axis.Axes')
            which_stat = display_option;
            display_option = ax;
            ax = gca;
        end
        
        display_option = validatestring(display_option, ...
            ["kstat" "corr" "MI" "shift"]);
        out = F.FitDiffs;
        
        switch display_option
            case "shift"
                % Show difference of means
                xx = abs(rad2deg(out.pdf_xcorr_angular_lag));
                yy = out.pdf_kstat_corr_aligned;
                F.gscatter_pat(ax, xx, yy, out.pat);
                title('PDF(M) v. PDF(D10)')
                xlabel('Shift (\circ)')
                ylabel('K-stat')
                xlim([-20 200]); xticks(-0:30:180); 
                ylim([0 1.1*max(yy)])
                
                % Summarize results
                disp(ax.Title.String);
                F.summarize_stat(ax.XLabel.String, xx, '', which_stat);
                F.summarize_stat(ax.YLabel.String, yy, '', which_stat);
                
            case "kstat"
                % Show difference of means
                xx = abs(rad2deg(out.Dmean));
                yy = out.pdf_kstat_corr_aligned;
                F.gscatter_pat(ax, xx, yy, out.pat);
                title('PDF(M) v. PDF(D10)')
                xlabel('Difference in means (\circ)')
                ylabel('K-stat')
                xlim([-20 200]); xticks(-0:30:180); 
                ylim([0 1.1*max(yy)])
                
                % Summarize results
                disp(ax.Title.String);
                F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang', which_stat);
                F.summarize_stat(ax.YLabel.String, yy, '', which_stat);
                
            case "MI"
                % Show difference of means
                xx = abs(rad2deg(out.Dmean));
                yy = out.mutual_info;
                F.gscatter_pat(ax, xx, yy, out.pat);
                title('PDF(M) v. PDF(D10)')
                xlabel('Difference in means (\circ)')
                ylabel('MI')
                xlim([-20 200]); xticks(-0:30:180); 
                ylim([0 1.1*max(yy)])
                
                % Summarize results
                disp(ax.Title.String);
                F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang', which_stat);
                F.summarize_stat(ax.YLabel.String, yy, '', which_stat);
                
            case "corr"
                % Show difference of means
                xx = abs(rad2deg(out.pdf_xcorr_angular_lag));
                yy = out.pdf_xcorr_coeff_normed;
                F.gscatter_pat(ax, xx, yy, out.pat);
                title('PDF(M) v. PDF(D10)')
                xlabel('Shift (\circ)')
                ylabel('peak CC (normalized)')
                xlim([-20 200]); xticks(-0:30:180); 
                ylim(quantile(yy, [0 1]) + .1*range(yy)*[-1 1])
                
                % Summarize results
                disp(ax.Title.String);
                F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang', which_stat);
                F.summarize_stat(ax.YLabel.String, yy, '', which_stat);
        
        end
        
        
        box off
        axis square
        

    end
    
    function ax = show_per_discharge_bias_vs_variance(F, ax, varargin)
        % ax = show_per_discharge_bias_vs_variance(F, ax=gca, ::method='median' ('min'/'max')::)
        % Shows patient gscatter of bias vs. variance
        if nargin < 2 || isempty(ax), ax = gca; end
        if ischar(ax), varargin = [ax varargin]; ax = gca; end
        
        method = 'Median';
        if mod(numel(varargin), 2)
            method = validatestring(varargin{1}, {'Min', 'Max', 'Median'}); 
        end
        switch method
            case 'Median'
                data_col = 1;
            case 'Max'
                data_col = 3;
            case 'Min'
                data_col = 2;
        end
        
        out = F.FitDiffs;
        
        xx = rad2deg(out.bias(:, data_col));
        yy = ((out.variance(:, data_col)));  
        F.gscatter_pat(ax, xx, yy, out.pat);
        title('Differences')
        xlabel([method ' bias (\circ)'])
        ylabel([method ' variance']);
        xlim([-20 200]); ylim([0 1]);
        xticks(-0:30:180); yticks(0:.25:1)
        box off
        axis square
        % Summarize results
        disp(ax.Title.String);
        F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang');
        F.summarize_stat(ax.YLabel.String, yy);
        
    end
    
    function ax = show_cc_between_methods(F, ax, use_mov_corr)
        % ax = show_cc_between_methods(F, ax=gca, use_mov_corr=false)
        % Shows the gscatter_pat of circular correlation (on y-axis; no
        % x-axis)
        
        if nargin < 2 || isempty(ax), ax = gca; end
        if nargin < 3 || isempty(use_mov_corr), use_mov_corr = false; end
        if islogical(ax); use_mov_corr = ax; ax = gca; end
        
        % Show CC between methods
        out = F.FitDiffs;
        switch use_mov_corr
            case true
                yy = out.CC(:, 1);
                F.gscatter_pat(ax, yy, out.pat);
                title(sprintf('Circular\ncorrelation'))
                % Summarize results
                disp(ax.Title.String)
                F.summarize_stat(ax.YLabel.String, yy);
                disp(out(yy < 0, ["pat" "sz" "CC"]))
                
            case false

                % Show CC between methods
                yy = out.corr_coeff;
                mask = out.corr_pval < .05;
                ln = F.gscatter_pat(ax, yy, out.pat);
                set(ln, 'markerfacecolor', 'none');
                hold on;
                F.gscatter_pat(ax, yy(mask), out.pat(mask));
                hold off
                
                % Summarize results
                disp(ax.Title.String)
                F.summarize_stat(ax.YLabel.String, yy);
                disp(out(yy < 0 & mask, ["pat" "sz" "corr_coeff" "corr_pval"]))
                
        end
        
        title(sprintf('Circular\ncorrelation'))
        ylabel('Corr coeff')
        ylim(quantile(yy, [0 1]) + .1*range(yy)*[-1 1])
        box off
        
        

    end
    
    function ax = show_max_dphi_dt(F, ax)
        % ax = show_max_dphi_dt(F, ax=gca)
        
        % Assumes only two methods.
        
        WIN = 10;
        OVERLAP = 9.9;
        step_ = WIN - OVERLAP;
        cellfun_ = @(varargin) cellfun(varargin{:}, 'uni', 0);
        
        if nargin < 2 || isempty(ax), ax = gca; end
        
        
        data = F.CompiledData;
        data = data(isfinite(data.dir), :);  % only use samples where phi is detected
        [G, pat, ~, mtc] = findgroups(data.patient, data.seizure, data.metric);

        % Compute dphi/dt
        dphi_dt = splitapply(@(dirs, times) ...
            {max(min(diff_phase(dirs) ./ diff(times), pi), -pi)}, ...  
            data.dir, data.time, G);
        sample_times = splitapply(@(t) {t(2:end)}, data.time, G);

        % Estimate mean dphi/dt (with ci) using WIN s windows
        test_times = cellfun_(@(t) t(1):step_:t(end), sample_times, 'uni', 0);
        for ii = 1:numel(test_times)
            temp = nan * test_times{ii};
            [~, locs] = min(abs(sample_times{ii}(:)' - test_times{ii}(:)));
            temp(locs) = dphi_dt{ii};
            dphi_dt{ii} = temp;
        end
        dphi_dt_mn = cellfun_(@(ddt, test_t) movmean(ddt, WIN, 'omitnan', 'samplepoints', test_t), dphi_dt, test_times);
        dphi_dt_std = cellfun_(@(ddt, test_t) movstd(ddt, WIN, 'omitnan', 'samplepoints', test_t), dphi_dt, test_times);
        N = cellfun_(@(ddt, test_t) movsum(isfinite(ddt), WIN, 'samplepoints', test_t), dphi_dt, test_times);
        dphi_dt_ci = cellfun_(@(std, N) 2*std./sqrt(N), dphi_dt_std, N);
        
        % Require 1 sample / second
        for ii = 1:numel(N)
            mask_ = N{ii} < WIN;
            dphi_dt_mn{ii}(mask_) = nan;
            dphi_dt_ci{ii}(mask_) = nan;
        end
        
%         % Estimate mean and bounds using circular stats
%         dphi_dt_mn = cellfun_(@(ddt, test_t, samp_t) ...  % this is the estimated mean dφ/dt at test_times in each patient
%             circ_mov_stat('mean', ddt, test_t, WIN, 'samplepoints', samp_t), ...
%             dphi_dt, test_times, sample_times, 'uni', 0); 
% 
%         dphi_dt_ci = cellfun_(@(ddt, test_t, samp_t) ...  % this is the ci around the estimated mean dφ/dt at test_times in each patient
%             circ_mov_stat('confmean', ddt, test_t, WIN, 'samplepoints', samp_t), ...
%             dphi_dt, test_times, sample_times, 'uni', 0); 
        
        % Convert to bounds
        bounds = cellfun_(@(mn, ci) mn(:) + [-ci(:) ci(:)], dphi_dt_mn, dphi_dt_ci);
        

        
%         % find points where ci do not include zero (i.e. dphi_dt
%         % significantly non-zero)
%         get_bounds = @(mn, ci) mn + [-ci ci];
%         inc_zero = @(bounds) isnan(bounds(:, 1)) | (bounds(:, 1) <= 0 & bounds(:, 2) >= 0);
%         is_significant = cellfun( ...
%             @(mn, ci) ~inc_zero(get_bounds(mn, ci)), ...
%             dphi_dt_mn, dphi_dt_ci, 'uni', 0);
        
%         % ... and set to 0 in both bounds and mn
%         fun = @(dat, is_sig) dat .* double(is_sig);
%         bounds = cellfun(@(mn, ci, is_sig) fun(get_bounds(mn, ci), is_sig), ...
%             dphi_dt_mn, dphi_dt_ci, is_significant, 'uni', 0);
%         dphi_dt_mn = cellfun(@(dat, is_sig) fun(dat, is_sig), ...
%             dphi_dt_mn, is_significant, 'uni', 0);
        
%         % Find max shift where ci are furthest from zero
%         get_ll_ = @(bounds) min(abs(bounds), [], 2);
%         peak_ll_ = @(ll) max(ll);
%         [ll_max, peak_non_zero] = cellfun_(@(bounds) peak_ll_(get_ll_(bounds)), bounds, 'uni', 0);
%         shift_max = cellfun_(@(mn, loc) mn(loc), dphi_dt_mn, peak_non_zero);
%         shift_max = abs(cat(1, shift_max{:}));
        
        % Find max dphi/dt
        shift_max = cellfun(@(x) max(abs(x)), dphi_dt_mn);
        
        isM1 = contains(F.Metrics, "M");  % Put M on the x-axis
        m1 = F.Metrics(isM1); m2 = F.Metrics(~isM1);
        assert(all(pat(mtc == m1) == pat(mtc == m2)));  % Make sure patient labels are the same for each metric
        
        % Make the plot
        xx = rad2deg(shift_max(mtc == m1));
        yy = rad2deg(shift_max(mtc == m2));
        plot([-20 200], [-20 200], '--', 'color', .5*[1 1 1]);
        hold on
        F.gscatter_pat(ax, xx, yy, pat(mtc == m1));
        hold off
        title('Peak d\phi/dt (\circ/s)')
        xlabel(m1)
        ylabel(m2)
        xlim([-20 200]); xticks(0:30:180); 
        ylim([-20 200]); yticks(0:30:180); 
        box off
        axis square
        % Summarize results
        disp(ax.Title.String);
        F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang');
        F.summarize_stat(ax.YLabel.String, deg2rad(yy), 'ang');
        
    end
    
        
    function ax = show_stat_with_range(F, ax, stat, units)
        % ax = show_stat_with_range(F, ax=gca, stat='bias')
        % Shows a stat for each seizure (e.g. bias, variance, cc) with its
        % range of observations from each 10s sliding window. Median is
        % shown with a marker; the line behind shows Q[0.05, 0.95])
        
        % Parse inputs
        if nargin < 2 || isempty(ax), ax = gca; end
        if nargin < 3 || isempty(stat), stat = 'bias'; end
        if nargin < 4 || isempty(units), units = ""; end
        if ischar(ax); units = stat; stat = ax; ax = gca; end
        sts = ax.NextPlot;
        
        % Choose a stat and get the dataset
        stat = validatestring(stat, {'Bias', 'Variance', 'CC', 'STD'});
        out = F.FitDiffs;
        
        out.bias = abs(out.bias); 
        out.std = out.mov_std_diffs;
        
        if units == "deg"
        
            out.bias = abs(rad2deg(out.bias));  % show bias in °
            out.std = rad2deg(circ_var2angle(out.mov_std_diffs));  % show std in °
            out.variance = rad2deg(circ_var2angle(2*out.variance));  % show variance in ° (convert to angular variance first)
            
        end
        
        
        % Make the plot
        ff = validatestring(stat, fields(out));
        yy = out.(ff)(:, 1);
        y_range = out.(ff)(:, 2:3);
        F.gscatter_pat(ax, yy, out.pat, 'range', y_range);
        switch stat
            case 'Bias'
                hold(ax, 'on')
                plot(ax, 0*yy, y_range(:, 2), 'go', 'markersize', 4);  % max bias
                plot(ax, 0*yy, yy, 'r.', 'markersize', 6);
                hold(ax, 'off')
                if units == "deg"
                    fprintf('%d seizures with max bias > 60°.\n', sum(abs(y_range(:, 2)) > 60));
                else
                    fprintf('%d seizures with max bias > π/3.\n', sum(abs(y_range(:, 2)) > pi/3));
                end
            case 'Variance'
                hold(ax, 'on')
                plot(ax, 0*yy, y_range(:, 1), 'bo', 'markersize', 4);  % min variance
                plot(ax, 0*yy, yy, 'r.', 'markersize', 6);
                hold(ax, 'off')
                if units == "deg"
                    fprintf('%d seizures with min var > 60°.\n', sum(abs(y_range(:, 1)) > 60));
                else
                    fprintf('%d seizures with min var > π/3.\n', sum(abs(y_range(:, 1)) > pi/3));
                end
            case 'STD'
                hold(ax, 'on')
                plot(ax, 0*yy, y_range(:, 1), 'bo', 'markersize', 4);  % min variance
                plot(ax, 0*yy, yy, 'r.', 'markersize', 6);
                hold(ax, 'off')
                if units == "deg"
                    fprintf('%d seizures with min std > 60°.\n', sum(abs(y_range(:, 1)) > 60));
                else
                    fprintf('%d seizures with min std > π/3.\n', sum(abs(y_range(:, 1)) > pi/3));
                end
        end
        
        % Prettify
        names = arrayfun(@(ll) string(ll.DisplayName), ax.Children);
        names(names == "") = [];
        xlim([0 12])
        xticks(ax, 1:11);  % there are 11 patients
        xticklabels(ax, flipud(names));
        ylabel(ax, stat)
        xlabel(ax, 'Patient');
        title('Median and Q[0.05, 0.95]')
        
        disp(stat)
        switch stat
            case 'Bias'
                if units == "deg"
                    ylim(ax, [-20 200]);
                    yticks(ax, 0:30:180)
                    ylabel(ax, [ax.YLabel.String ' (°)'])
                else
                    ylim(ax, deg2rad([-20 200]));
                    yticks(ax, deg2rad(0:30:180))
                    ylabel(ax, ax.YLabel.String)
                    polar_label(ax, 'y');
                end
                
                F.summarize_stat([stat ' (med)'], deg2rad(yy), 'ang');
                F.summarize_stat([stat ' (min)'], deg2rad(y_range(:, 1)), 'ang');
                F.summarize_stat([stat ' (max)'], deg2rad(y_range(:, 2)), 'ang');
                
            case 'Variance'
                if units == "deg"
                    ylim(ax, [-20 200]);
                    yticks(ax, 0:30:180)
                    ylabel(ax, [ax.YLabel.String ' (°)'])
                    
                    F.summarize_stat([stat ' (med)'], deg2rad(yy), 'ang');
                    F.summarize_stat([stat ' (min)'], deg2rad(y_range(:, 1)), 'ang');
                    F.summarize_stat([stat ' (max)'], deg2rad(y_range(:, 2)), 'ang');
                else
                    ylim(ax, [0 1]);
                    ylabel(ax, ax.YLabel.String)
                    
                    F.summarize_stat([stat ' (med)'], (yy), '');
                    F.summarize_stat([stat ' (min)'], (y_range(:, 1)), '');
                    F.summarize_stat([stat ' (max)'], (y_range(:, 2)), '');
                end
                
                
                
            case 'STD'
                if units == "deg"
                    ylim(ax, [-20 200]);
                    yticks(ax, 0:30:180)
                    ylabel(ax, [ax.YLabel.String ' (°)'])
                else
                    ylim(ax, [0 1]);
                    ylabel(ax, ax.YLabel.String)
                end
                
                F.summarize_stat([stat ' (med)'], deg2rad(yy), 'ang');
                F.summarize_stat([stat ' (min)'], deg2rad(y_range(:, 1)), 'ang');
                F.summarize_stat([stat ' (max)'], deg2rad(y_range(:, 2)), 'ang');
                
            otherwise
                F.summarize_stat([stat ' (med)'], yy);
                F.summarize_stat([stat ' (min)'], y_range(:, 1));
                F.summarize_stat([stat ' (max)'], y_range(:, 2));
        end
        
        ax.NextPlot = sts;
        
            
    end
    
    function allP_pdf_distance2(F)
        

        Nax = 16;  % 14 is ideal for a row of plots
%         r = 3; c = ceil(Nax/r);
        h = figure('name', 'pdf_distance', 'units', 'inches', ...
            'position', [0 0 Nax+1 3.25] * .67);  
        
        Ts = tiledlayout(h, 1, Nax, 'tilespacing', 'none');
        
        
        % Show comparison of distributions
        ax = nexttile(Ts, [1 3]);
        F.show_distribution_comparison(ax, "shift");
        
        % Show range of bias and variance for each seizure
        aB = nexttile(Ts, [1 5]);
        F.show_stat_with_range(aB, 'bias', 'deg');
        
        aV = nexttile(Ts, [1 5]);
        F.show_stat_with_range(aV, 'var');

        set([aB aV], 'xlabel', [], 'xticklabelrotation', 45)
        
        
        % Show max shift
        ax = nexttile(Ts, [1 3]);
        F.show_max_shift(ax);
        
        ln = findobj(h, 'type', 'line');
        set(ln, 'markersize', 4);
        
%         % Show CC between methods
%         ax = nexttile(Ts, [1 2]);
%         F.show_cc_between_methods(ax, false);

%         % Show per discharge comparison (distribution of differences)
%         % Bias and variance (median and max) on sliding 10 s windows
%         ax = nexttile(Ts, [1 3]);
%         ax = F.show_per_discharge_bias_vs_variance(ax, 'median');
%         
%         ax = nexttile(Ts, [1 3]);
%         ax = F.show_per_discharge_bias_vs_variance(ax, 'max');        
       
        
        
%         % Show Kstat of differences from VM(mu_{d\phi}, kappa_{d\phi})
%         % i.e. how VM are the differences?
%         ax = nexttile(Ts, [1 2]);
%         yy = out.mov_kstat;
%         F.gscatter_pat(ax, yy, out.pat);
%         title('Diffs v. VM')
%         ylabel('Avg. Kuiper stat')
%         ylim(quantile(yy, [0 1]) + .1*range(yy)*[-1 1])
%         box off
%         % Summarize results
%         summarize(ax.Title.String, yy);
%         
%         
%         
%         % Show variance (peak and median)
%         ax = nexttile(Ts, [1 3]);
%         [~, ~, aa] = F.allP_directionality;
%         copyobj(aa.Children, ax);
%         set(ax, 'xlabel', aa.XLabel, 'ylabel', aa.YLabel)
%         close(aa.Parent);
%         lines = findobj(ax, 'type', 'line');
%         % Convert to variance (1 - DI)
%         for ll = lines'
%             ll.XData = 1 - ll.XData;
%             ll.YData = 1 - ll.YData;
%         end
%         ln_medians = arrayfun(@(ln) strcmpi(ln.MarkerFaceColor, 'none'), lines);
%         xx_meds = cat(2, lines(ln_medians).XData);
%         yy_meds = cat(2, lines(ln_medians).YData);
%         xx_pks = cat(2, lines(~ln_medians).XData);
%         yy_pks = cat(2, lines(~ln_medians).YData);
%         title(ax, sprintf('Var(%s) v. Var(%s)', F.Metrics(1), F.Metrics(2)))
%         % Summarize results
%         disp(ax.Title.String)
%         disp('medians')
%         summarize(ax.XLabel.String, xx_meds);
%         summarize(ax.YLabel.String, rmmissing(yy_meds));
%         disp('peaks')
%         summarize(ax.XLabel.String, xx_pks);
%         summarize(ax.YLabel.String, rmmissing(yy_pks));
%         legend(ax, chillins);
        
        
        
%         % Show directionality index
%         ax = nexttile(Ts, [1 2]);
%         mask1 = contains(mtc, "D"); % assuming comparison of an M and D metric and I want the D on the x-axis
%         mask2 = contains(mtc, "M"); 
%         DI = splitapply(@nanmean, data.dir_ind_05, G);
%         xx = DI(mask1);
%         yy = DI(mask2);
%         assert(all(pat(mask1) == pat(mask2)));
%         pat = pat(mask1); sz = sz(mask1);
%         F.gscatter_pat(ax, xx, yy, pat);
%         xlabel(unique(mtc(mask1)));
%         ylabel(unique(mtc(mask2)));
%         title('DI')
%         % Summarize results
%         disp(ax.Title.String)
%         summarize(ax.XLabel.String, xx);
%         summarize(ax.YLabel.String, yy);
%         mask = yy > xx;
%         fprintf('DI_M > DI_D in:\n')
%         cellfun(@disp, compose('%s %d', pat(mask), sz(mask)));
%        

%         
        
        
        
        
        
        
%         % Show Kstat of differences from VM(mu_{d\phi}, kappa_{d\phi})
%         ax = nexttile(Ts, [1 2]);
%         yy = out.mov_kstat;
%         F.gscatter_pat(ax, yy, out.pat);
%         title('Diffs v. VM')
%         ylabel('Kuiper stat')
%         ylim(quantile(yy, [0 1]) + .1*range(yy)*[-1 1])
%         box off
%         % Summarize results
%         disp(ax.Title.String)
%         summarize(ax.YLabel.String, yy);
        
        

        
%         % Show distribution difference of means
%         ax = nexttile(Ts, [1 3]);
%         xx = rad2deg(out.pdf_xcorr_angular_lag);
%         yy = rad2deg(out.Dmean);
%         F.gscatter_pat(ax, xx, yy, out.pat);
%         title('Overall shift')
%         xlabel('CC shift (\circ)')
%         ylabel('Difference in means (\circ)')
%         xlim([-180 180]); xticks(-180:60:180); 
%         ylim([-180 180]); yticks(-180:60:180); 
%         box off
%         axis square
%         % Summarize results
%         disp(ax.Title.String);
%         summarize_ang(ax.XLabel.String, deg2rad(xx));
%         summarize_ang(ax.YLabel.String, deg2rad(yy));
%         
%         
%         
%         % Show distribution distance (correlation distance and Kuiper stat)
%         ax = nexttile(Ts, [1 3]);
%         xx = out.pdf_xcorr_coeff;
%         yy = 1 - out.pdf_kstat;
%         F.gscatter_pat(ax, xx, yy, out.pat);
%         xlim([0 1])
%         ylim([0 1])
%         title('Distribution similarity');
%         ylabel("1 - Kstat");
%         xlabel('Cross correlation');
%         axis square
%         box off
%         % Summarize results
%         disp(ax.Title.String);
%         summarize(ax.XLabel.String, xx);
%         summarize(ax.YLabel.String, yy);
%         
        
        
        
        
        % Only show one legend
        lgd = findobj(gcf, 'type', 'legend');
        set(lgd(2:end), 'visible', 'off');
        lgd(1).Location = 'eastoutside';


        F.print(h, F.prefix_better(''));
    end

    
    function allP_pdf_distance(F)
        
        out = F.FitDiffs;
        out.Dmean_diffs = rad2deg(out.Dmean_diffs);
        out.Dmode_d = rad2deg(out.pdf_xcorr_angular_lag);


        fields = {'Dmean_diffs', 'Dstd_diffs', 'Dmode_d', 'mutual_info'};
        labels = {'Per discharge difference (\circ)', ...
            'Per discharge std', ...
            'Difference in modes (\circ)', ... 
            'MI'};
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
                    yline(.3, 'r--', 'linewidth', 1)
            end
        end

        lgd = findobj(gcf, 'type', 'legend');
        set(lgd(2:end), 'visible', 'off');
        lgd(1).Location = 'eastoutside';

%         set(Ts.Children, 'fontsize', 11)

        F.print(h, F.prefix_better(''));
    end

    function allP_main_iw(F, min_electrodes)
        % allP_main_iw(F, min_electrodes=50)
        if nargin < 2 || isempty(min_electrodes), min_electrodes = 50; end
        
        iw_all = BVNY.get_iw_main;
        Nch = structfun(@(s) numfinite(s.template.template(s.main_wave, :)), iw_all);
        mask = Nch >= min_electrodes;
        h = figure('name', 'iw_templates'); fullwidth(true);
        T = tiledlayout(h, 4, 8);
        
        fields = fieldnames(iw_all);
        [px, py] = ndgrid(1:10, 1:10);
        pos = [px(:) py(:)];
        for ff = string(fields(mask)')
            ax = nexttile(T);
            tpl = squeeze( ...  % TOA
                iw_all.(ff).template.template(iw_all.(ff).main_wave, :, :) ...
                );
            M = MaxDescent;
            M.Position = pos;
            M.t0 = 0;
            M.Data = tpl;
            M.Name = ff;
            M.plot2D(ax);
            
            % Convert TOA to seconds (if you do this before now, the
            % speed in the legends will be off)
            sc = findobj(ax, 'type', 'scatter');
            sc.CData = sc.CData / 1e3;  % Convert from ms -> s
            ax.CLim = ax.CLim / 1e3;
            cb = ax.Colorbar;
            title(cb, 'TOA [s]');
            
            
            axis square
            title(strrep(ff, '_', ' '))
            xlim([0 11])
            ylim([0 11])
            xticks([])
            yticks([])
            box on
        end
        
        % Set colomap to gray (light to dark)
        % This is mirroring the recoloring |F.F.show_corr_dir_mismatches|.
        % Hoping to easily distinguish between IW and method-specific TW
        gray_ = .5*[1 1 1];
        len = length(ax.Colormap);
        cmap = make_diverging_colormap( ...
            [brighten(gray_, .9); brighten(gray_, -.9)], gray_, len);
        set(findobj(T, 'type', 'axes'), 'colormap', cmap);
        
        % Make quivers colored
        qq = findobj(T, 'type', 'quiver');
        set(qq, 'color', 1-lines(1));
        
        % Put legends outside
        lgd = findobj(T, 'type', 'legend');
        set(lgd, 'location', 'southoutside');
        
        % Set all scatter dots to the same size
        sc = findobj(T, 'type', 'scatter');
        set(sc, 'sizedata', M.get_dot_size(ax));
        
        
        F.print(h, F.prefix_better(''));
        
    end
    
    function ZZallP_main_iw(F, min_electrodes)  % ZZ'd 6/1/21
        % allP_main_iw(F, min_electrodes=50)
        if nargin < 2 || isempty(min_electrodes), min_electrodes = 50; end
        
        iw_all = BVNY.get_iw_main;
        Nch = structfun(@(s) numfinite(s.template.template(s.main_wave, :)), iw_all);
        mask = Nch >= min_electrodes;
        h = figure('name', 'iw_templates'); fullwidth(true);
        T = tiledlayout(h, 'flow');
        
        fields = fieldnames(iw_all);
        for ff = string(fields(mask)')
            tpl = squeeze(iw_all.(ff).template.template(iw_all.(ff).main_wave, :, :));  % TOA
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
    
    function allP_dir_v_rho(F, ax, metric, min_electrodes)
        % Correlation v. angle (relative to IW angle)
        if nargin < 2 || isempty(ax), ax = []; end
        if nargin < 3 || isempty(metric), metric = "M"; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 50; end
        if ~isaxes(ax) && nargin == 3, min_electrodes = metric; metric = ax; ax = []; end
        if ~isaxes(ax) && nargin == 2, metric = ax; ax = []; end
        
        if isempty(ax)
            h = figure('name', 'rho v phi', 'position', [0 0 4 4]);
            ax = axes(h);
        else
            h = ax.Parent;
            while ~isa(h, 'matlab.ui.Figure')  % ... in case the axis is in a tiledlayout
                h = h.Parent;
            end
        end
        
        data = F.CompiledData;
        mask = strcmpi(data.metric, metric) ...
            & isfinite(data.dir) ...
            & data.nchannels >= min_electrodes; % ...
        data = data(mask, :);
        
        % Some renaming
        data.dir0 = data.dir;
        data.dir_d = rad2deg(fix_angle(data.dir - data.iw_angle));
        
        % Find outliers (i.e. points that are corr+ and d180 or corr- and
        % d0
        [~, ~, resid] = regress(data.rho, [ones(size(data.rho)) abs(data.dir_d)]);
        big_r = abs(normalize(resid)) > 2;
        
        ln_bigR = F.gscatter_pat(ax, ...
            data.dir_d(big_r), data.rho(big_r), data.patient(big_r));
        hold(ax, 'on');
        ln_all = F.gscatter_pat(ax, data.dir_d, data.rho, data.patient);
        hold(ax, 'on');
        xlabel(ax, '\phi [\circ]')
        ylabel(ax, '\rho')
        lgd = legend(ax, 'location', 'eastoutside');
        title(ax, sprintf("Direction v. corr (%s)", metric))
        xlim(ax, [-180 180])
        ylim(ax, [-1 1])
        xticks(ax, -180:60:180)
        xline(ax, 0); yline(ax, 0)
        lgd.String(contains(lgd.String, 'data')) = '';
        axis(ax, 'square')
%         ln = findobj(ax, 'type', 'line');
        set(ln_all, 'markersize', 1);
        set(ln_bigR, 'markersize', 4);
        
        F.print(h, F.prefix_better(metric));
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
            M = F.get_wave_fit(F.get_file(ss));

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
            M = F.get_wave_fit(F.get_file(ss));

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
        M0.direction_scatter(ax, 'markerfacecolor', [0 0 0], 'sizedata', 2);
        hold(ax, 'on');
        Mnew.direction_scatter(ax, 'markerfacecolor', [1 0 0], 'sizedata', 2);
        hold(ax, 'off');
        
        sc = get(ax, 'children');
        
    end
    function dirs = get.SlidingDirs(F)
        % Computes the estimated mean direction in sliding MEAN_WIN=10s windows
        % (OVERLAP=9.9s).
        % Returns a structure with the estimated means with CI for each TOA
        % method in F.Metrics.
        % Takes a little bit to compute and used to show the max shift and
        % to count the large shifts to better to compute and store.
        
        % If stored, get value and return
        if ~isempty(F.SlidingDirs), dirs = F.SlidingDirs; return; end
        
        MEAN_WIN = 10;
        SHIFT_WIN = 10;
        OVERLAP = 9.9;
        STEP_ = MEAN_WIN - OVERLAP;
        
        
        
        data = F.CompiledData;
        data = data(isfinite(data.dir), :);  % only use samples where phi is detected
        
        % Compute everything groupwise (grouped by seizure and metric)
        [G, pat, sz, mtc] = findgroups(data.patient, data.seizure, data.metric);
        
        wng = warning;
        warning('off', 'circ_confmean:requirementsNotMet');
        
        dirs = splitapply(@(x) {x}, data.dir, G);
        sample_times = splitapply(@(t) {t}, data.time, G);
        test_times = cellfun(@(t) (t(1):STEP_:t(end))', sample_times, 'uni', 0);
%         test_times = sample_times;
        dir_mn = cellfun(@(dirs, test_t, samp_t) ...  % get moving circular mean on 10 s windows
            (circ_mov_stat('mean', dirs, test_t, MEAN_WIN, 'samplepoints', samp_t)), ...
            dirs, test_times, sample_times, 'uni', 0);
        dir_ci = cellfun(@(dirs, test_t, samp_t) ...  % get moving circular mean on 10 s windows
            (circ_mov_stat('confmean', dirs, test_t, MEAN_WIN, 'samplepoints', samp_t)), ...
            dirs, test_times, sample_times, 'uni', 0);
        N = cellfun( ...  % number of sample points in the interval
            @(samp_t, test_t) sum(abs(samp_t(:) - test_t(:)') <= MEAN_WIN/2), ...
            sample_times, test_times, 'uni', 0);
        for ii = 1:numel(N)
            mask_ = N{ii} < MEAN_WIN/2;
            dir_mn{ii}(mask_) = nan;
            dir_ci{ii}(mask_) = nan;
        end
        
        % Get the max shift occuring within a SHIFT_WIN s interval centered at
        % each time
        filter_times = @(tt) double(abs(tt - tt') <= SHIFT_WIN/2);  % Creates a diagonal matrix with ones indicating nearby times 
        angle_diff = @(mn) abs(fix_angle(mn - mn'));  % Creates a diagonal matrix of pairwise angular distances
        shifts = cellfun( @(mn, t) ...
            max( ...  
                angle_diff(mn) .* filter_times(t), ...  
            [], 2), dir_mn, test_times, 'uni', 0);
        
        
        % I think return a table of everything similar to CompiledData
        % Convert all table variables to vectors
        Gind = arrayfun(@(ii) repmat(ii, numel(test_times{ii}), 1), ...
            1:length(test_times), 'uni', 0);
        Gind = cat(1, Gind{:});
        patient = pat(Gind);
        seizure = sz(Gind);
        metric = mtc(Gind);
        time = cat(1, test_times{:});
        mean = cat(1, dir_mn{:});
        ci = cat(1, dir_ci{:});
        intervalN = cat(2, N{:})';
        shifts = cat(1, shifts{:});
        
        % Put it all together ]
        dirs = table(patient, seizure, metric, time, mean, ci, intervalN, shifts);
        
        dirs.Properties.UserData = struct(...
            'MeanWin', MEAN_WIN, 'ShiftWin', SHIFT_WIN, 'Overlap', OVERLAP);
        F.SlidingDirs = dirs;
        
                
        warning(wng)
        
    end
    
    
    function ax = show_max_shift(F, ax)
        % ax = show_max_shift(F, ax=gca)
        % Shows the largest shift (in degrees) detected by each method over
        % any 10 second window. 
        
        
        WIN = 10;
        fprintf('Computing max shift occuring within a %d second interval.\n', WIN);
        
        if nargin < 2 || isempty(ax), ax = gca; end
                
        data = F.SlidingDirs;
        [G, pat, ~, mtc] = findgroups(data.patient, data.seizure, data.metric);
        
        
        shift_max = splitapply(@max, data.shifts, G);
        
        
        isM1 = contains(F.Metrics, "M");  % Put M on the x-axis
        m1 = F.Metrics(isM1); m2 = F.Metrics(~isM1);
        assert(all(pat(mtc == m1) == pat(mtc == m2)));  % Make sure patient labels are the same for each metric
        
        % Make the plot
        xx = double(mtc == m1) + 1;
        yy = rad2deg(shift_max);

        x2 = [xx(mtc == m1), xx(mtc == m2)];
        y2 = [yy(mtc == m1), yy(mtc == m2)];
        pp = [pat(mtc == m1), pat(mtc == m2)];

        % show lines between methods
        for ii = 1:size(x2, 1)
            F.gscatter_pat(ax, x2(ii, :)', y2(ii, :)', pp(ii, :)');
            hold(ax, 'on')
        end
        hold(ax, 'off');
        ln = findobj(ax, 'type', 'line');
        set(ln, 'linestyle', '-', 'linewidth', .5);
        names = arrayfun(@(ll) string(ll.DisplayName), ln);
        ln = ln(~ismember(names, ["data1" "data2" ""]));
        names = names(~ismember(names, ["data1" "data2" ""]));
        [nn, uu] = unique(names);
        [pp, pu] = unique(F.SeizureInfo.patient);
        [~, so] = sort(pu);
        map = containers.Map(string(pp), so);
        legend(ln(uu(arrayfun(@(p) map(p), nn))));

        title('Max shift')
        xlabel('')
        xlim([0.5 2.5]);
        xticklabels([m1 m2]);
        xticks([1 2]);
        ylabel('Magnitude (°)');
        ylim([-20 200]); yticks(0:30:180); 
        box off
        axis square
        % Summarize results
        disp(ax.Title.String);
        F.summarize_stat(m1, deg2rad(y2(:, 1)), 'ang');
        F.summarize_stat(m2, deg2rad(y2(:, 2)), 'ang');
        
       
    end
    
    function ax = ZZshow_max_shift(F, ax)  % ZZ'd 5/29/21
        % ax = show_max_shift(F, ax=gca)
        % Shows the largest shift (in degrees) detected by each method over
        % any 10 second window. 
        % First computes circular moving mean over 10s windows, then finds
        % the max range (again, sliding 10s windows)
        % Assumes only two methods.
        
        WIN = 10;
        OVERLAP = 9.9;
        STEP_ = WIN - OVERLAP;
        
        if nargin < 2 || isempty(ax), ax = gca; end
        
        
        data = F.CompiledData;
        data = data(isfinite(data.dir), :);  % only use samples where phi is detected
        [G, pat, ~, mtc] = findgroups(data.patient, data.seizure, data.metric);
        
        dirs = splitapply(@(x) {x}, data.dir, G);
        sample_times = splitapply(@(t) {t}, data.time, G);
        test_times = cellfun(@(t) t(1):STEP_:t(end), sample_times, 'uni', 0);
%         test_times = sample_times;
        dir_mn = cellfun(@(dirs, test_t, samp_t) ...
            (circ_mov_stat('mean', dirs, test_t, WIN, 'samplepoints', samp_t)), ... % get moving circular mean on 10 s windows
            dirs, test_times, sample_times, 'uni', 0);
        dir_ci = cellfun(@(dirs, test_t, samp_t) ...
            (circ_mov_stat('confmean', dirs, test_t, WIN, 'samplepoints', samp_t)), ... % get moving circular mean on 10 s windows
            dirs, test_times, sample_times, 'uni', 0);
        N = cellfun( ...
            @(samp_t, test_t) sum(abs(samp_t(:) - test_t(:)') <= WIN/2), ...  % number of sample points in the interval
            sample_times, test_times, 'uni', 0);
        for ii = 1:numel(N)
            mask_ = N{ii} < WIN/2;
            dir_mn{ii}(mask_) = nan;
            dir_ci{ii}(mask_) = nan;
        end
        
        shift_max = cellfun(@(x, times) ...
            max(abs(angle(exp(1j * (x - x')))) ...  % Get the largest angular shift
                .* double(abs(times - times') <= WIN/2), ...  % ... within WIN/2 s
            [], 'all'), ...
            dir_mn, test_times);
        
        isM1 = contains(F.Metrics, "M");  % Put M on the x-axis
        m1 = F.Metrics(isM1); m2 = F.Metrics(~isM1);
        assert(all(pat(mtc == m1) == pat(mtc == m2)));  % Make sure patient labels are the same for each metric
        
        % Make the plot
        xx = rad2deg(shift_max(mtc == m1));
        yy = rad2deg(shift_max(mtc == m2));
        plot([-20 200], [-20 200], '--', 'color', .5*[1 1 1]);
        hold on
        F.gscatter_pat(ax, xx, yy, pat(mtc == m1));
        hold off
        title('Max shift (\circ)')
        xlabel(m1)
        ylabel(m2)
        xlim([-20 200]); xticks(0:30:180); 
        ylim([-20 200]); yticks(0:30:180); 
        box off
        axis square
        % Summarize results
        disp(ax.Title.String);
        F.summarize_stat(ax.XLabel.String, deg2rad(xx), 'ang');
        F.summarize_stat(ax.YLabel.String, deg2rad(yy), 'ang');
        
    end

    
    
    function h = allS_dir_v_time(F, min_electrodes)
        % Shows the direction v time for each seizure using each method
        % sc = allS_dir_v_time(F, min_electrodes=50)
        if nargin < 3 || isempty(min_electrodes), min_electrodes = 50; end
        

        % Get the data
        data = F.CompiledData;
        mask = isfinite(data.dir) ...  % ... finite direction 
            & data.nchannels >= min_electrodes;  % ... and # electrodes recruited to IW
        data = data(mask, :);
        
        % Center all times to the IW
        data.dir_time_offset = data.time - data.iw_center;
        
        % Put units in degrees and center to IW
        data.dir = fix_angle(data.dir - data.iw_angle);  % center to IW
        data.dir_deg = rad2deg(data.dir);

        % Group by patient
        [G, pat, sz] = findgroups(data.patient, data.seizure);
        pat_num = findgroups(pat);
        sz_num = cell2mat(splitapply(@(x) {findgroups(x)}, sz, pat_num));
        c = max(pat_num); r = max(sz_num);
        
        % Create the figure and tiled layout
        h = figure('name', 'dir_v_time', 'position', [0 0 2*c 1.25*r]);
        T = tiledlayout(h, r, c, 'tilespacing', 'compact', 'padding', 'none');
        lw = .5;
%         gray_ = .85*[1 1 1];

        % For each seizure create a scatter plot of the rho value for each
        % discharge time. Show all rho in gray; highlight significant rho
        % in color; use different colors for each seizure
        tile_num = sub2ind(fliplr(T.GridSize), pat_num, sz_num);

        for ii = 1:max(G)
            
            % Show all discharges colored by metric
            
            ax = nexttile(T, tile_num(ii));
            xx = data.dir_time_offset(G == ii);
            yy = data.dir_deg(G == ii);
            for mm = F.Metrics
                mask_ = data.metric(G == ii) == mm;
                plot(ax, xx(mask_), yy(mask_), '.', 'markersize', 4, ...
                    'color', F.Style.(mm).color, 'displayname', mm);
                hold(ax, 'on');
            end
            hold(ax, 'off')
            
            
            % Show the seizure onset time for each seizure
            sz_onset = -unique(data.iw_center(G == ii));
            xline(ax, sz_onset, ':', 'linewidth', lw, 'displayname', 'Sz onset');
            
            % prettify
            yticks([-180 -90 0 90 180]);
            yticklabels(["-180°" "" "0°" "" "180°"])
            ylim([-200 200])
            title(sprintf('%s %d', pat(ii), sz(ii)))
            xline(ax, 0, 'linewidth', lw, 'displayname', 'IW passage');  % highlight IW time
            hold off
            grid on
            ylabel('');
            xlabel('Time [s]');
%             xticklabels('');
            
        end
        
        % prettify
        xlabel('Time [s]')
        ax = T.Children;
        linkaxes(ax, 'xy')
        set(ax, 'xtickmode', 'auto');
        ax(1).XTickLabelMode = 'auto';
        
        % Move axes from tiledlayout straight to figure so you can put a
        % legend in the bottom corner
        ax = copyobj(ax, h);
        delete(T);
        lgd = legend(ax(end), 'location', 'south');
        lgd.Position(2) = 0;
        
        
        F.print(h, F.prefix_better(sprintf('%s_%s_%d', F.Metrics, min_electrodes)));
    end
    
    function sc = allP_dir_v_time(F, metric, min_electrodes)
        % Shows the direction v time for each patient (different seizures
        % identified by color)
        % sc = allP_corr_v_time(F, metric='M', thresh=5e-2, min_electrodes=0)
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(min_electrodes), min_electrodes = 40; end
        

        % Get the data
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) ...  % mask by metric,
            & isfinite(data.dir) ...  % ... finite direction 
            & data.nchannels >= min_electrodes;  % ... and # electrodes recruited to IW
        data = data(mask, :);
        
        % Center all times to the IW
        data.dir_time_offset = data.time - data.iw_center;
        
        % Put units in degrees
        data.dir_deg = rad2deg(data.dir);
%         data.dir_deg = rad2deg(angle(exp(1j*(data.dir - data.iw_angle))));

        % Group by patient
        [G, pat] = findgroups(data.patient);
        r = numel(unique(pat));
        
        % Create the figure and tiled layout
        h = figure('name', 'dir_v_time', 'position', [0 0 3 r]);
        T = tiledlayout(h, r, 1);
        iw_yline = @(yval, color) yline(yval, 'color', color);
        onset_xline = @(xval, color) xline(xval, 'color', color);
        lw = .5;
%         gray_ = .5*[1 1 1];

        % For each patient create a scatter plot of the rho value for each
        % discharge time. Show all rho in gray; highlight significant rho
        % in color; use different colors for each seizure
        for ii = 1:numel(pat)
            
            % Show all discharges colored by seizure number
            ax = nexttile(T, ii);
            xx = data.dir_time_offset(G == ii);
            yy = data.dir_deg(G == ii);
            cc = data.seizure(G == ii);  % color by seizure
            sc = gscatter(ax, xx, yy, cc, [], [], 6, 'off');

            % Get the color corresponding to each seizure
            cols = flipud(arrayfun(@(l) l.Color, sc, 'uni', 0));
            
            % Show the IW direction for each seizure
            iw_angles = rad2deg(unique(data.iw_angle(G == ii)));
            ln = arrayfun( ...
                @(ii) iw_yline(iw_angles(ii), cols{ii}), 1:numel(cols));
            set(ln, 'linewidth', lw);  % don't know why, but putting this in the iw_yline function doesn't work
            
            % Show the seizure onset time for each seizure
            sz_onset = -unique(data.iw_center(G == ii));
            ln = arrayfun( ...
                @(ii) onset_xline(sz_onset(ii), cols{ii}), 1:numel(cols));
            set(ln, 'linewidth', lw);
            
            % prettify
            yticks([-180 0 180]);
            title(sprintf('%s', pat(ii)))
            xline(0)
            hold off
            grid on
            ylabel('\phi');
            xlabel('');
            xticklabels('');
            
        end
        
        % prettify
        xlabel('Time [s]')
        ax = T.Children;
        linkaxes(ax, 'xy')
        set(ax, 'xtickmode', 'auto');
        ax(1).XTickLabelMode = 'auto';
        F.print(h, F.prefix_better(metric));

    end

    
    function sc = allP_corr_v_time(F, metric, thresh, min_electrodes)
        % Shows the correlation v time for each patient (different seizures
        % identified by color)
        % sc = allP_corr_v_time(F, metric='M', thresh=5e-2, min_electrodes=0)
        if nargin < 2 || isempty(metric), metric = 'M'; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 40; end
        
        % Create ksdensity local fun
        XI = linspace(-1, 1, 100);
        ksdensity_ = @(xx) ksdensity(xx, XI, ...
            'bandwidth', .1, 'boundarycorrection', 'log', ...
            'support', [-1 1]);

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
        h = figure('name', 'corr_v_time', 'position', [0 0 4.5 r]);
        T = tiledlayout(h, r, 4);
        

        % For each patient create a scatter plot of the rho value for each
        % discharge time. Use different colors for each seizure
        for ii = 1:numel(pat)
            
            % Show all discharges in gray
            ax_sc = nexttile(T, [1 3]); 
            set(ax_sc, 'tag', 'scatter', 'nextplot', 'replacechildren');
            xx = data.rho_time_offset(G == ii);
            yy = data.rho(G == ii);
            cc = findgroups(data.seizure(G == ii));  % color by seizure
            sc = scatter(ax_sc, xx, yy, 1, cc, 'filled');
            ax_sc.Colormap = hsv(max(cc));

            % highlight significant rho in color
            hold on
            mask = data.rho_pval(G == ii) < thresh;
            sc = scatter(ax_sc, xx(mask), yy(mask), 4, cc(mask), 'filled');
            
            % prettify
            ylim([-1 1]);
            title(sprintf('%s', pat(ii)))
            xline(0)
            hold off
            grid on
            ylabel('\rho');
            
            % Get PDF of rho values
            ax_hist = nexttile(T, [1 1]); 
            set(ax_hist, 'tag', 'hist', 'nextplot', 'replacechildren');
            pdf = splitapply(@(rho) ksdensity_(rho), yy, cc);
            plot(ax_hist, XI, pdf');
            ylabel(ax_hist, 'PDF');
            set(ax_hist, 'colororder', hsv(max(cc)))
            
            
        end
        
        % prettify
        xlabel(ax_sc, 'Time [s]')
        xlabel(ax_hist, '\rho');
        set(findobj(h, 'tag', 'hist'), 'yaxislocation', 'right');
        linkaxes(findobj(h, 'tag', 'scatter'), 'x')        
        F.print(h, F.prefix_better(metric));

    end
    
    function dataR = allP_hist2d_tw_v_iw_effectsize(F, metric, thresh, min_electrodes, rho_or_dir)
        % dataR = allP_hist2d_tw_v_iw_effectsize(F, ...
        %       metric="M", thresh=5e-2, min_electrodes=50, rho_or_dir='dir')
        % 
        % Shows an estimate of the probability that a discharge will be a
        % 0°|IW (180°|IW) discharge through time.
        % 
        % For t in [seizure_onset:STEP=2:seizure_termination], resample
        % from the set of discharges in a 4 second (2*HALFWIN) window
        % centered at t and observed in seizures where an IW is detected on
        % at least <min_electrodes>.

        
        if nargin < 2 || isempty(metric), metric = "M"; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 50; end
        if nargin < 5 || isempty(rho_or_dir), rho_or_dir = 'dir'; end
        rho_or_dir = validatestring(rho_or_dir, ["rho" "dir"]);
        metric = string(metric);

        Nouter = 1e4;
        HALFWIN = 2;
        STEP = 2;
        ptest = .05; % Show CI to this level
        
        if rho_or_dir == "rho"
            dnames = ["Corr" "Anticorr"];
        else
            dnames = ["\phi < 90" "\phi > 90"];
        end
        
        h = figure('name', 'hist2d_dir', 'position', [0 0 2 1.5] * 2);
        data = F.CompiledData;
        mask = strcmp(data.metric, metric) ...
            & data.nchannels > min_electrodes;
        data = data(mask, :);
        
        
        % Use direction instead of correlation
        if rho_or_dir == "rho"
            phi = sign(data.rho) .* double(data.rho_pval < thresh);
            sample_times = data.rho_time_offset;
        else  % rho_or_dir == "dir"
            % remove non-finite directions (i.e. where direction could not be estimated)
%             data = data(isfinite(data.dir), :);  
            dir = angle(exp(1j*(data.dir - data.iw_angle)));
            phi = double(abs(dir) <= pi/2) - double(abs(dir) > pi/2);
            sample_times = data.time - data.iw_center;
        end
        
        % Group data by patient; get indices;
        % determine the lowest number of discharges in each
        % patient. Use this as the resample number
        G = findgroups(data.patient);
        N = min(splitapply(@numel, G, G));
        
        trange = quantile(sample_times, [0 1]);
        tbins = linspace(trange(1), trange(2), diff(trange)/STEP + 1);  % Divide times into ~2s intervals

        
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
        T = tiledlayout(h, 3, 1);
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
        for jj = 1:2  % Show corr and anti corr proportions

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
                'facecolor', gray_, 'linestyle', 'none', 'facealpha', .5); %#ok<NASGU>

            % Show the mean
            plot(ax, tbinsI, mnI(:, jj), 'color', color, ...
                'linewidth', 2, 'displayname', dnames(jj))
            hold(ax, 'off');

            % highlight IW time 
            xline(ax, 0)

            % Prettify
            grid(ax, 'on')
            axis(ax, 'tight');
            title(ax, dnames(jj));

        end
        
        
        %%% Compare the corr/anticorr distributions
        childs = get(findobj(T, 'type', 'axes'), 'children');
        ax = nexttile(T, 3);
        
        % Shade regions where distributions differ
        xx = tbinsI';
        mask = QhiI(:, 2) <= QloI(:, 1) | QhiI(:, 1) <= QloI(:, 2);
        ylo = min(QhiI, [], 2); yhi = max(QloI, [], 2);
        ylo(~mask) = yhi(~mask); %yhi(~mask) = mnI(~mask, 1);
%         ylo(~mask) = 0; yhi(~mask) = 0;
        yy = [ylo; flipud(yhi)];
        xx = [xx; flipud(xx)];
        xx(isnan(yy)) = [];
        yy(isnan(yy)) = [];
        
        pp = patch(ax, xx, yy, 1, 'facecolor', [0 0 0], ...  % Using the gray face color is too hard to see
            'facealpha', .8, 'linestyle', 'none'); %#ok<NASGU>
        
        % Copy the objects in the first two axes into the third axis
        cc = cellfun(@(cc) copyobj(cc, ax), childs, 'uni', 0); %#ok<NASGU>
        delete(findobj(ax, 'type', 'area'));
        patches = findobj(ax, 'type', 'patch');  % these are the CI
        lns = findobj(ax, 'type', 'line');  % these are the medians
        cl = findobj(ax, 'type', 'constantline');  % delete one of the constant lines
        delete(cl(1));
        
        % Darken the patch in back
        darker_ = brighten(color, -.75);  % Darken the color
        set(patches(2), 'facecolor', darker_);
        set(lns(2), 'color', darker_);
        ax.XLim = T.Children(2).XLim;
        
        
        % Add a title
        title(ax, 'Comparison');

        
%         if rho_or_dir == "rho", ylim(ax, [0 1]); else, ylim(ax, [0 .5]); end
%         xlim(ax, [-30 50]);  % Nothing significant beyond this
        xlabel(ax, 'Time [s]');
        ylabel(T, 'Proportion');
        
        ax = findobj(T, 'type', 'axes');
        linkaxes(ax, 'xy');
        set(ax, 'box', 'on', 'ylim', [0 1]);
        lgd = legend(ax(2), {'95%CI', 'CI>0', 'median'}, 'location', 'eastoutside'); %#ok<NASGU>

        % Summary statement
        fields = ["pos" "neg"];
        for jj = 1:2
            [~, loc] = max(QloI(:, jj));
            rho_at_t.(fields(jj)) = mnI(loc, jj);
            peak_t.(fields(jj)) = tbinsI(loc);
            rho_ci_at_t.(fields(jj)) = [QloI(loc, jj) QhiI(loc, jj)];
        end
        
        fprintf(['Using this method, we find that the highest ' ...
            'probability of seeing a correlated discharge occurs at ' ...
            't=%0.2f s (p_corr = %0.2f, [%0.2f, %0.2f]; median, 95%%CI); ' ...
            'the highest probability of seeing an anti-correlated ' ...
            'discharge occurs at t=%0.2f s (p_anti = %0.2f, [%0.2f, %0.2f]; '...
            'median, 95%%CI; Figure 4C)\n'], ...
            peak_t.pos, rho_at_t.pos, rho_ci_at_t.pos, ...
            peak_t.neg, rho_at_t.neg, rho_ci_at_t.neg)
    

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
    
    function [yy, xx, pat] = allP_percent_sigest_rho_v_dir(F, metric, thresh, min_electrodes)
        % [yy, xx, pat] = allP_percent_sigest_rho_v_dir(F, ...
        %       metric="M", thresh=5e-2, min_electrodes=50 )
        % Plots the percent of the discharges that have strong
        % relationships with the main IW template. i.e. Percentage of
        % correlations with p-value less than thresh for each patient and
        % seizure.
        if nargin < 2 || isempty(metric), metric = "M"; end
        if nargin < 3 || isempty(thresh), thresh = 5e-2; end
        if nargin < 4 || isempty(min_electrodes), min_electrodes = 50; end
        
        metric = validatestring(metric, F.Metrics);
        
        figure('units', 'inches', 'position', [0 0 2 1.25] * 2, ...
            'name', 'TW-IW corr');
        
        data = F.CompiledData;
        mask = ...  % only use data where
            ismember(data.metric, metric) ...  % the TOA method matches
            & isfinite(data.rho_pval) ...  % the correlation could be computed (dir will also be nan here)
            & data.nchannels >= min_electrodes;  % the IW was detected on at least min_electrodes
        data = data(mask, :);

        [G, pat] = findgroups(data.patient, data.seizure);
        

        xx = splitapply(@(x) mean(isfinite(x))*100, data.dir, G);
        lx = strcat(metric, " [% plane != (0,0)]");

        
        yy = splitapply(@(x) mean(x < thresh)*100, data.rho_pval, G);
        F.gscatter_pat(xx, yy, pat);
        ax = gca;
        ylabel(strcat(metric, " [% corr != 0]"))
        xlabel(lx)
        title({'TW detected'; ...
            sprintf('(p<5e%0.0f)', log10(thresh/5))});
        
        lgd = legend('location', 'eastoutside');
        axis equal
        xlim([50 100]);
        ylim([0 100]);
        
        % Add a dashed line to show y=x
        line(ax, [0 100], [0 100], 'linestyle', '--', 'color', .5*[1 1 1]);
        lgd.String(end) = []; % remove line from legend
        
        F.print(F.prefix_better(metric));
        
        % Summary statement
        fprintf(['%s:\nIn %d/%d seizures, at least 50%% ', ...
            'of discharges have a non-zero correlation with '...
            'the IW pattern (median=%0.1f%%, range=[%0.1f%%,%0.1f%%]).\n'], ...
            metric, sum(yy >= 50), numel(yy), quantile(yy, [.5 0 1]));
        
        fprintf(['This is in contrast to the plane fitting ' ...
            'method, where we find significant fits ' ...
            'less often (median=%0.1f%%, range=[%0.1f%%,%0.1f%%]).\n'], ...
            quantile(xx, [.5 0 1]));
        
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
        ff = 'rho_pval';  % dir
        if m2 == ""
%             yy = splitapply(@(x) sum(x < thresh), data.rho_pval(mask(m1)), G(mask(m1)));
            yy = splitapply(@numel, data.(ff)(mask(m1)), G(mask(m1)));
            ly = strcat(m1, " [#]");
            tag = m1;
        else
            yy = splitapply(@(x) mean(x < thresh)*100, data.(ff)(mask(m2)), G(mask(m2)));
            ly = strcat(m2, " [% corr\neq0]");
            tag = sprintf('%s_%s', m1, m2);
        end
        xx = splitapply(@(x) mean(x < thresh)*100, data.(ff)(mask(m1)), G(mask(m1)));
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
        
        % Summary statement
        fprintf(['In %d/%d seizures, at least 50%% ', ...
            'of discharges have a non-zero correlation with '...
            'the IW pattern (median=%0.1f%%, range=[%0.1f%%,%0.1f%%]).\n'], ...
            sum(xx >= 50), numel(xx), quantile(xx, [.5 0 1]));
        
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
    
    function [xx, yy, ax] = allP_directionality(F, metrics, win)
        % Shows the mean directionality index over win s intervals
        
        data = F.CompiledData;
        if nargin < 2, metrics = F.Metrics; end
        if nargin < 3, win = 10; end
        
        figure('units', 'inches', 'position', [0 0 3 2]);
        mask = contains(data.metric, metrics) ...
            & isfinite(data.dir);
        data = data(mask, :);

        m2 = metrics{1}; m1 = metrics{2};

        G = findgroups(data.patient, data.seizure, data.metric);
        data.num_finite_dirs = nan(size(data.time));
        for ii = 1:max(G)
            mask = G == ii;
            times = data.time(mask);
            N = arrayfun(@(tt) sum( abs(times - tt) < win/2 ), times);  % Use a 5 second window
            data.num_finite_dirs(mask) = N;
        end
        field = sprintf('dir_ind_%02d', win);
        data.(field)(data.num_finite_dirs < win) = nan;  % only look where we get at least 10 samples
        
        [G, pat, sz] = findgroups(data.patient, data.seizure);

        mask =@(mm) strcmpi(data.metric, mm);

        % show min and median
        xx = splitapply(@(x) quantile(x, [0 .5]), data.dir_ind_10(mask(m1)), G(mask(m1)));
        yy = splitapply(@(x) quantile(x, [0 .5]), data.dir_ind_10(mask(m2)), G(mask(m2)));
        
        
        sc = F.gscatter_pat(xx(:, 1), yy(:, 1), pat);
        hold on
        sc = F.gscatter_pat(xx(:, 2), yy(:, 2), pat);
        set(sc, 'markerfacecolor', 'none', 'markersize', 4);
        hold off
                
%         disagree = F.Disagree;
%         dd = cellfun(@(x) find(strcmpi(x, compose('%s %d', pat, sz))), disagree);
%         
%         ax = gca;
%         hold(ax, 'on')
%         plot(ax, xx(dd, :), yy(dd, :), 'k*', 'markersize', 6);
%         hold(ax, 'off')
        
        title('Directionality index');
        xlabel(m1)
        ylabel(m2)
        xlim([0 1]);
        ylim([0 1]);
        axis square
        legend off
        ax = gca;

        F.print(F.prefix_better(sprintf('%s_%s_w%d', m1, m2, win)));
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
    
    function allP_iw_nchannels(F)
        data = F.get_iw_table(0);  % load the iw stats
%         mask = data.wave_num == data.main_wave;
%         data = data(mask, :);
        
        data.disagree = ismember( ...
            compose('%s %d', string(data.patient), data.seizure), ...
            F.Disagree.Row ...
            );
        
        % Describe results
        [G, pat, sz] = findgroups(string(data.patient), data.seizure);
        Niwc = splitapply(@numel, G, G);
        res = table(pat, sz, Niwc);
        disp(res)
        summary(res)
        fprintf('%d seizures with at least 1 IWc\n', max(G));
        fprintf('%d seizures with at least 2 IWc\n', sum(Niwc > 1));
        
        
        
        h = figure('name', 'iw_nchannels', 'position', [0 0 .6*4 1.5] * 2);
        T = tiledlayout(h, 1, 4);
        
        
        ax = nexttile(T, [1 2]);
        % show all waves in outline
        ln = F.gscatter_pat(data.nchannels, data.patient);
        set(ln, 'markerfacecolor', 'none');
        hold(ax, 'on')
        for ll = ln', ll.ZData = -1*ones(size(ll.XData)); end
        


        % show disagrees with asterisks
        mask = data.disagree & data.wave_num == data.main_wave;
        ln = F.gscatter_pat(data.nchannels(mask), data.patient(mask));
        set(ln, 'marker', '*', 'color', [0 0 0], 'markersize', 6);
        for ll = ln', ll.ZData = 1*ones(size(ll.XData)); end

        
        % show main wave in color
        [G, pat] = findgroups(data.patient, data.seizure);
        yy = splitapply(@(a, b, c) c(a == b), data.wave_num, data.main_wave, data.nchannels, G);
        ln = F.gscatter_pat(yy, pat);             
        for ll = ln', ll.ZData = 0*ones(size(ll.XData)); end

        hold off
        set(findobj(T, 'type', 'line'), 'linewidth', 1)
        title('Nchan')
        ylabel('[#]')
        legend('location', 'westout')
        
        
        ax = nexttile(T, [1 2]);
        dat = data(data.wave_num == data.main_wave, :);
        dsub = data(data.wave_num ~= data.main_wave, :);
        max_bin = ceil(max(data.nchannels)/10) * 10;
        histogram(ax, dsub.nchannels, 10:10:max_bin)
        hold(ax, 'on')
        histogram(ax, dat.nchannels, 10:10:max_bin)
        hold(ax, 'off');
        set(ax.Children, 'edgecolor', [1 1 1], 'edgealpha', .5);
        xlabel(ax, 'Nchan'); ylabel(ax, 'IW candidates [counts]')
        legend(ax, 'Secondary', 'Main')
        title(ax, 'Channels recruited to IW')
        
        F.print(h, F.prefix_better(''));
    end
    
    function [yy, ax] = allP_iw_stats(F, min_electrodes)
        % Show info about detected IW for each patient and seizure
        
        if nargin < 2, min_electrodes = 0; end
        PCUTOFF = -8;
        SPEEDCUTOFF = 8;
        
        fields = ["fr_peak_mu", "iw_fwhm_mu", "prop_sp", "p"];        
        h = figure('name', 'iw_stats', 'position', [0 0 .6*(numel(fields)+1) 1.5] * 2);

        data = F.get_iw_table(0);  % load the iw stats
        data.prop_sp(data.p >= .05) = nan;  % only show speeds where the fit was significant;
        data.p = log10(data.p/5);  % rescale p-value so that it's in terms of p=5eXX
        data.width = data.iw_fwhm_mu .* data.prop_sp;
        
        % rescale very low pvalues
        lo_p = data.p < PCUTOFF;
        data.p(lo_p) = rescale(data.p(lo_p), PCUTOFF-.5, PCUTOFF-.05);
        
        % rescale very high speeds
        hi_speed = data.prop_sp > SPEEDCUTOFF;
        data.prop_sp(hi_speed) = rescale(data.prop_sp(hi_speed), SPEEDCUTOFF + .1, SPEEDCUTOFF + 1);

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
            'iw_fwhm_mu', {{'Duration', '[s]', 'linear'}}, ...
            'width', {{'Width', '[mm]', 'linear'}}, ...
            'crossing_time', {{'Crossing time', '[s]', 'linear'}});

        data.disagree = ismember( ...
            compose('%s %d', string(data.patient), data.seizure), ...
            F.Disagree.Row);

        T = tiledlayout(gcf, 1, numel(fields));
        for ff = fields
            ax = nexttile(T);
            dat = data;
            
            % Separate mains and secondaries
            mains = data.wave_num == data.main_wave;
            xx = double(mains & data.nchannels >= min_electrodes);  % x=0 -> secondary; x=1 -> main
            yy = dat.(ff);
            % Summarize mean with STD
            F.summarize_stat(ff, rmmissing(yy(xx == 1)), [], 'std');
            F.ttest2_(rmmissing(yy(xx == 0)), rmmissing(yy(xx == 1)), ...
                'vartype', 'equal');
            

            % show secondary waves in outline
            ln = F.gscatter_pat(xx(~mains), yy(~mains), dat.patient(~mains));
            set(ln, 'markerfacecolor', 'none');
            hold(ax, 'on')
            for ll = ln', ll.ZData = -1*ones(size(ll.XData)); end


            % show disagrees with asterisks
            mask = dat.disagree & dat.wave_num == dat.main_wave;
            ln = F.gscatter_pat(xx, yy(mask), dat.patient(mask));
            set(ln, 'marker', '*', 'color', [0 0 0], 'markersize', 6);
            for ll = ln', ll.ZData = 1*ones(size(ll.XData)); end

            % show main wave in color
            [G, pat] = findgroups(dat.patient, dat.seizure);
            yy = splitapply(@(a, b, c) c(a == b), dat.wave_num, dat.main_wave, dat.(ff), G);
            ln = F.gscatter_pat(xx(mains), yy, pat);             
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
            xlim(quantile(xx, [0 1]) + .5*[-1 1])
            xticks([0 1])
            xticklabels(["Exc" "Inc"])
        end
        lgd = findobj(T, 'type', 'legend');
        set(lgd, 'visible', 'off');
        set(lgd(1), 'visible', 'on', 'location', 'eastoutside')
        lgd(1).String(contains(lgd(1).String, 'data')) = [];
        F.print(h, F.prefix_better(''));
    end
    
    function allP_compare_metrics(F, metrics)
        % Print all to a giant figure
        if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
        

        N = numel(F.SeizureInfo.patient);
        r = 5; c = ceil(N/r);
        
        h = figure('name', 'allP_compare_metrics', ...
            'position', [0 0 c r] .* [0 0 4.5 2], 'visible', 'off');
        h.Units = 'Normalized';
        
        % Get position coordinates (divide the grid)
        [bottom, left] = ndgrid(1-1/r:-1/r:0, 0:1/c:1-1/c);
        w = 1/c; ht = 1/r;
        
        % add margins
        left = left + .1*w;
        bottom = bottom + .1*ht;
        
        subplot_pos = @(ii) [left(ii) bottom(ii) .8*w .6*ht];
        
        for ii = 1:N
            pat = F.get_fname(ii);
            htemp = F.compare_metrics2(pat, metrics, false);
            lgd = findobj(htemp, 'type', 'legend');
            if ii == 1, vis = 'on'; else, vis = 'off'; end
            set(lgd, 'visible', vis, 'location', 'northeast');
            T = copyobj(htemp.Children, h);
            set(T, 'position', subplot_pos(ii));
            close(htemp)
        end
        
        F.print(h, F.prefix_better(sprintf('%s_%s', metrics)));
        close(h);
    end

    function allDisagree_compare_metrics(F, metrics)
        if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
        

        disagree = F.Disagree;
        N = numel(disagree);
        r = 5; c = ceil(N/r);
        
        h = figure('name', 'disagrees_compare_fits', ...
            'position', [0 0 c r] .* [0 0 2 1] * 2);
        h.Units = 'Normalized';
        
        % Get position coordinates (divide the grid)
        [bottom, left] = ndgrid(1-1/r:-1/r:0, 0:1/c:1-1/c);
        w = 1/c; ht = 1/r;
        
        % add margins
        left = left + .1*w;
        bottom = bottom + .1*ht;
        
        subplot_pos = @(ii) [left(ii) bottom(ii) .8*w .6*ht];
        
        for ii = 1:N
            pp = disagree{ii};
            pat = strrep(pp, ' ', '_Seizure');
            htemp = F.compare_metrics(pat, metrics, false);
            lgd = findobj(htemp, 'type', 'legend');
            if ii == 1, vis = 'on'; else, vis = 'off'; end
            set(lgd, 'visible', vis, 'location', 'northeast');
            T = copyobj(htemp.Children, h);
            set(T, 'position', subplot_pos(ii));
            close(htemp)
        end
        
        F.print(F.prefix_better(sprintf('%s_%s', metrics)));
    end
    
    function h = compare_metrics2(F, pat, metrics, print_flag)
        % Plot direction rasters from metrics for patient in pat
        
        
        if nargin < 2 || isempty(pat), pat = "c7_Seizure1"; end
        if isnumeric(pat), pat = string(F.SeizureInfo.name{pat}); end
        if nargin < 3, metrics = string(F.Metrics(1:2)); end
        if nargin < 4, print_flag = true; end
        
        pat = string(validatestring(pat, F.SeizureInfo.name));
        sz_num = find(strcmp(pat, F.SeizureInfo.name));
        
        h = figure('name', 'compare_fits', 'position', [0 0 3*1.5+1 2]);
        T = tiledlayout(h, 2, 3, 'TileSpacing', 'compact');
        
        if pat == "c7_Seizure1", rotate_by = -pi/2; 
        elseif contains(pat, "MG49"), rotate_by = pi;
        else, rotate_by = 0; 
        end
            
        S1 = [1 2];  % times
        S2 = [2 1];  % pdfs
        
        m1 = metrics(1);
        m2 = metrics(2);
        
        fit1 = F.get_wave_fit(F.get_fname(sz_num));
        fit2 = fit1.(m2);
        fit1 = fit1.(m1);
        
        % Get the diffs
        if contains(m2, "D")
            fitD = copy(fit1);
            fitD.Direction = fit1.diff(fit2);
        else
            fitD = copy(fit2);
            fitD.Direction = fit2.diff(fit1);
        end
        fitD.RotateBy = 0;
        
        t1 = fit1.time;
        t2 = fit2.time;
        tD = fitD.time;
        
        fit1.RotateBy = rotate_by;
        fit2.RotateBy = rotate_by;
        
        d1 = fit1.Direction;
        d2 = fit2.Direction;
        dD = fitD.Direction;
        
        
        
        % Direction v. time plots
        colors = arrayfun(@(mm) F.Style.(mm).color, [m1 m2], 'uni', 0);
        ax1 = nexttile(T, S1);
        fit1.plot_dirs_with_ci(ax1, 'color', colors{1}, 'displayname', m1, 'linewidth', 1);
        hold(ax1, 'on')
        fit2.plot_dirs_with_ci(ax1, 'color', colors{2}, 'displayname', m2, 'linewidth', 1);
        scatter(ax1, t1, rad2deg(d1), 2, ones(size(t1)), 'filled');
        scatter(ax1, t2, rad2deg(d2), 2, 2 * ones(size(t2)), 'filled')
        hold(ax1, 'off')
        set(ax1, 'colormap', brighten(cat(1, colors{:}), -.5), ...
            'xlabel', [], 'xticklabel', []);
        ylabel('Dirs')
%         xlabel(ax1, 'Time [s]');
        
        % Show the differences
        ax2 = nexttile(T, S1);
        fitD.plot_dirs_with_ci(ax2, 'linewidth', 1);
        hold(ax2, 'on')
        scatter(ax2, tD, rad2deg(dD), 2, ones(size(tD)), ...
            'filled', 'displayname', 'Obs');
        pp = findobj(ax2, 'type', 'patch');
        ax2.Colormap = pp.FaceColor;
        hold(ax2, 'off')
        ln = findobj(ax2, 'type', 'line');
        pp = findobj(ax2, 'type', 'patch');
        scD = findobj(ax2, 'type', 'scatter');
%         legend(ax2, [ln(1) pp(1) sc(1)], 'location', 'eastoutside');
        ylabel(ax2, 'Diffs');
        xlabel(ax2, 'Time [s]');       
        
        aa = [ax1 ax2];
        linkaxes(aa, 'x');
        set(aa, 'ytick', [-180 0 180], 'yaxislocation', 'left', ...
            'yticklabel', ["-180°" "0°" "180°"], 'title', [])
        sc = findobj(aa, 'type', 'scatter');
        set(sc, 'sizedata', 2);
        
        % PDF plot
        ax = nexttile(T, S2);
        axD = ax;
%         axD = nexttile(T, S2);
        [f1, xi] = circ_ksdens(d1);
        f2 = circ_ksdens(d2, xi);
        fd = circ_ksdens(dD, xi);
        
        if 0  % Get sliding window distributions?
            
            % Get fit1 sliding window distributions
            all_f1 = nan(numel(t1), numel(xi)); %#ok<UNRCH>
            for ii = 1:numel(t1)
                mask = abs(t1 - t1(ii)) < WIN/2;
                if numfinite(d1(mask)) < WIN, continue; end
                all_f1(ii, :) = circ_ksdens(d1(mask), xi);
            end

            % Get fit2 sliding window distributions
            all_f2 = nan(numel(t2), numel(xi));
            for ii = 1:numel(t2)
                mask = abs(t2 - t2(ii)) < WIN/2;
                if numfinite(d2(mask)) < WIN, continue; end
                all_f2(ii, :) = circ_ksdens(d2(mask), xi);
            end

            % Get diffs sliding window distributions
            all_diffs = nan(numel(diff_times), numel(xi));
            for ii = 1:numel(diff_times)
                mask = abs(diff_times - diff_times(ii)) < WIN/2;
                if numfinite(diffs(mask)) < WIN, continue; end
                all_diffs(ii, :) = circ_ksdens(diffs(mask), xi);
            end
        
            if 0  % Show quantile patches?
                getQ = @(dat) quantile(dat, [.025 .975]);
                q1 = getQ(all_f1);
                q2 = getQ(all_f2);
                qd = getQ(all_diffs);
                p1 = patch(ax, rad2deg([xi(:); flipud(xi(:))]), [q1(1, :)'; flipud(q1(2, :)')], 1, ...
                    'facecolor', F.Style.(m1).color, 'facealpha', .5, 'linestyle', 'none');
                p2 = patch(ax, rad2deg([xi(:); flipud(xi(:))]), [q2(1, :)'; flipud(q2(2, :)')], 1, ...
                    'facecolor', F.Style.(m2).color, 'facealpha', .5, 'linestyle', 'none');
                p3 = patch(axD, rad2deg([xi(:); flipud(xi(:))]), [qd(1, :)'; flipud(qd(2, :)')], 1, ...
                    'facecolor', .15 * [1 1 1], 'facealpha', .3, 'linestyle', 'none');
            end
        
            f1 = nanmedian(all_f1);
            f2 = nanmedian(all_f2);
            fd = nanmedian(all_diffs);
        end
        
        hold([ax, axD], 'on'); 
        lw = 1;
        p3 = plot(axD, rad2deg(xi), fd, ':', 'color', .15*[1 1 1], ...
            'displayname', 'Diffs', 'linewidth', lw); 
        p1 = plot(ax, rad2deg(xi), f1, 'displayname', m1, ...
            'color', brighten(F.Style.(m1).color, -.5), 'linewidth', lw);
        p2 = plot(ax, rad2deg(xi), f2, 'displayname', m2, ...
            'color', brighten(F.Style.(m2).color, -.5), 'linewidth', lw);
        hold([ax axD], 'off');
        legend(ax, [p1 p2 p3 ln(1) pp(1) scD(1)], 'location', 'eastoutside');
        xticks(ax, -180:90:180);
        ax.XTickLabel = cellfun(@(tk) sprintf('%s°', tk), ax.XTickLabel, 'uni', 0);
        xlabel(axD, 'Direction')
        title(ax, 'PDF');
        ylim(ax, [0 1.1*max([f1, f2, fd])]);
%         xticklabels(ax, []);
        linkaxes([ax axD], 'x');
        grid(ax, 'on');
        
        
        % Labels
        title(ax1, F.SeizureInfo.display_names{sz_num}, 'fontsize', 11);
        
        if print_flag

            % Print
            pat = F.SeizureInfo.patient{sz_num};
            tag = sprintf('%s_%s_%s', pat, m1, m2);
            F.print(h, F.prefix_better(tag));
            
        end
        
    end

    function h = compare_metrics(F, pat, metrics, print_flag)
        % Plot direction rasters from metrics for patient in pat
        if nargin < 2 || isempty(pat), pat = "c7_Seizure1"; end
        if isnumeric(pat), pat = string(F.SeizureInfo.name{pat}); end
        if nargin < 3, metrics = string(F.Metrics(1:2)); end
        if nargin < 4, print_flag = true; end
        
        sz_num = find(strcmpi(pat, F.SeizureInfo.name));
        
        h = figure('name', 'compare_fits', 'position', [0 0 2.5 1] * 2);
        T = tiledlayout(h, 3, 6, 'TileSpacing', 'compact');
        
        if pat == "c7_Seizure1", rotate_by = -pi/2; 
        elseif contains(pat, "MG49"), rotate_by = pi;
        else, rotate_by = 0; 
        end
            
        S1 = [1 3];
        S2 = [3 3];
        
        m1 = metrics(1);
        m2 = metrics(2);
        
        fit1 = F.get_wave_fit(F.get_fname(sz_num));
        fit2 = fit1.(m2);
        fit1 = fit1.(m1);
        
        fit1.RotateBy = rotate_by;
        fit2.RotateBy = rotate_by;
        
        % Get the diffs
        if numel(fit1.time) < numel(fit2.time)
            diffs = fit1.diff(fit2);
            diff_raster = @(ax) fit1.direction_raster(ax, diffs);
        else
            diffs = fit2.diff(fit1);
            diff_raster = @(ax) fit2.direction_raster(ax, diffs);
        end
        
        % Density plot
        ax = nexttile(T, 4, S2);
        [f1, xi] = circ_ksdens(fit1.Direction);
        f2 = circ_ksdens(fit2.Direction);
        f3 = circ_ksdens(diffs, xi);
        
        p3 = plot(ax, rad2deg(xi), f3, 'k:', 'displayname', 'Diffs'); 
        hold(ax, 'on'); 
        p1 = plot(ax, rad2deg(xi), f1, 'displayname', m1, 'color', F.Style.(m1).color);
        p2 = plot(ax, rad2deg(xi), f2, 'displayname', m2, 'color', F.Style.(m2).color);
        hold(ax, 'off');
        legend(ax, [p1 p2 p3], 'location', 'northeastoutside');
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
        diff_raster(ax3);
        ylabel('Diff');
        xlabel('Time [s]');
        title(ax3, '');
        
        aa = [ax1 ax2 ax3];
        linkaxes(aa, 'x');
        set(aa, 'ytick', [-180 0 180], 'yaxislocation', 'left')
        sc = findobj(aa, 'type', 'scatter');
        set(sc, 'sizedata', 2);
        
        % Labels
        title(T, F.SeizureInfo.display_names{sz_num}, 'fontsize', 11);
        
        if print_flag

            % Print
            pat = F.SeizureInfo.patient{sz_num};
            tag = sprintf('%s_%s_%s', pat, m1, m2);
            F.print(h, F.prefix_better(tag));
            
            % Make a colorwheel
            pax = fit1.colorwheel;
            F.print(pax.Parent, F.prefix_better('colorwheel'));
            close(pax.Parent);
        end
        
    end
    
    function data = compile_fits(F, interval, dir_index_win)
        % int is the interval surrounding the IW to assign phase 1 (early)
        % or phase 3 (late)
        
        if nargin < 2 || isempty(interval), interval = Inf; end  
        if nargin < 3 || isempty(dir_index_win), dir_index_win = 5; end
        USE_MID = true;  % compute the intervals from the midpoint of the IW (or from the bounds)
        
        % Load the IW and TW data
        W = F.WaveFits;
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
                
                
%                 [dir, sp] = fit.discharge_directions;  % This was an
%                 attempt at identifying isolated discharges (not relevant
%                 now that M is tested at unique discharges)
                dir = fit.Direction;
                sp = fit.Magnitude;
                data_temp = table( ...
                    dir, ...
                    fit.time, ...
                    sp, ...
                    [nan; fit.dtheta], ...
                    fit.directionality_index(1), ...
                    fit.directionality_index(2), ...
                    fit.directionality_index(5), ...
                    fit.directionality_index(10), ...
                    'VariableNames', {'dir', 'time', 'speed', ...
                    'dtheta', ...
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
        F.allP_pdf_distance2;
        F.allP_proportion_finite;
        
        % Make examples for one agree patient and one disagree patient
        for pp = F.ExampleSeizures, F.compare_metrics2(pp); end
        
        % only show stats for seizures where IW was detected on at least 30
        % electrodes. There are some cases where an IW is detected on
        % fewer electrodes but this is less likely to actually be an IW
        MIN_ELEC = 50;
        F.allP_iw_stats(MIN_ELEC); 
        

        
        %%% F3: TW DIRECTION WRT IW DIRECTION %%%
        % Here we show that there is no clear evidence to support the IW
        % hypothesis, but also not enough against it to reject
        
        % make the combined direction v. time with IW plots for each
        % patient
        F.allS_dir_v_time;
        F.dir_dist_by_phase();  % no longer used, but still in Illustrator so keep updating
        F.make_polarhist;  % show the pre-/post-IW polar histograms for each patient and seizure
        
        % Make the stick plots and polar histograms for all patients
        % inf is the interval surrounding the IW; I chose 40
        % (min_electrodes) since the main wave of most seizures appears on
        % at least 50 (40/50 makes no difference)
        F.pre_post_direction_stick_plots([], MIN_ELEC);  
        F.pre_post_direction_hist([], MIN_ELEC);
        
        
        %%% F4: TW/IW CORRELATION %%%
        % Here we show that [...TBD]
        
        % Show correlation v. time (relative to IW)
        F.allP_corr_v_time('M', 5e-2, MIN_ELEC);
        
        % Show that the correlation measure is giving us more information
        % and how it relates to the direction measure
        F.allP_percent_sigest_rho_v_dir('M', 5e-2, MIN_ELEC);
        F.allP_dir_v_rho("M");
        
        % Test for times of significance (this one takes a while)
        F.allP_hist2d_tw_v_iw_effectsize('M', 5e-2, MIN_ELEC, 'rho');
        
        % Show how many TW have a strong relationship to the main IW. Only
        % show these numbers for IW that appear on more than 50 electrodes.
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
    
    function [hyp, p, ci, stats] = ttest2_(x, y, varargin)
        % Wrapper for ttest2. Prints some information to stdout
        [hyp, p, ci, stats] = ttest2(x, y, varargin{:});
        ci_fun = @(x) 1.96*std(x)/sqrt(numel(x))*[-1 1] + mean(x);
        mn = @mean;
        fprintf('Group 1 mean, [95%%CI]: %0.3g, [%0.3g, %0.3g]\n', mn(x), ci_fun(x));
        fprintf('Group 2 mean, [95%%CI]: %0.3g, [%0.3g, %0.3g]\n', mn(y), ci_fun(y));
        fprintf('Reject null: %s, p=%0.3g\n', string(logical(hyp)), p);
        fprintf('Estimated difference in means: [%0.3g, %0.3g]\n', ci);
    end
   
    function summarize_stat(str, xx, method, stat)
        % methods are '', 'ang', or 'halfang'; stat is 'ci', 'std', 'range'
        % '': normal
        % 'ang': circular mean
        % 'halfang': half normal converted to °
        % Using stat='ci' means estimating the distribution mean with CI;
        % using stat='std' means estimating the sample mean with STD; using
        % stat='range' returns the median and range of the observations
        % stat='iqr' returns the median and IQR
        
        if nargin < 3 || isempty(method), method = ''; end
        if nargin < 4 || isempty(stat), stat = "iqr"; end
        method = validatestring(method, {'', 'ang', 'halfang'});
        stat = validatestring(stat, ["ci", "std", "range", "iqr"]);
        str = sprintf('%8s', str);  % Align messages
        
        switch method
            case ''
                if stat == "ci"
                    fprintf('%s: mean, [95%%CI] = %0.3g, [%0.3g, %0.3g]\n', ...
                        str, mean(xx), mean(xx) + 2/sqrt(numel(xx))*std(xx) * [-1 1]);
                elseif stat == "std"
                    fprintf('%s: mean (STD) = %0.3g (%0.3g)\n', ...
                        str, mean(xx), std(xx));
                elseif stat == "range"
                    fprintf('%s: median [range] = %0.3g [%0.3g, %0.3g]\n', ...
                        str, median(xx), quantile(xx, [0 1]));
                elseif stat == "iqr"
                    fprintf('%s: median [IQR] = %0.3g [%0.3g]\n', ...
                        str, nanmedian(xx), iqr(xx));
                end
            case 'ang'
                if stat == "ci"
                    fprintf('%s: circ mean, [95%%CI] = %0.3f°, [%0.3f°, %0.3f°]\n', ...
                        str, rad2deg(circ_mean(rmmissing(xx))), ...
                        rad2deg(circ_mean(rmmissing(xx)) + circ_confmean(rmmissing(xx)) * [-1 1]));
                elseif stat == "std"
                    fprintf('%s: circ mean (circ var) = %0.3f° (%0.3f)\n', ...
                        str, rad2deg(circ_mean(rmmissing(xx))), ...
                        circ_var(rmmissing(xx)));
                elseif stat == "range"
                    md = circ_median(xx);
                    Q = quantile(fix_angle(xx - md), [0 1]);
                    Q = fix_angle(Q + md);
                    fprintf('%s: median [range] = %0.3g° [%0.3g°, %0.3g°]\n', ...
                        str, rad2deg(md), Q);
                elseif stat == "iqr"
                    md = circ_median(rmmissing(xx));
                    Q = quantile(fix_angle(xx - md), [0.25 .75]);
                    Q = fix_angle(Q + md);
                    fprintf('%s: median [IQR] = %0.3g° [%0.3g°]\n', ...
                        str, rad2deg(md), rad2deg(fix_angle(diff(Q))));
                end
            case 'halfang'
                % ... you're never going to use this... not updating
                sig2mean = @(sigma) sigma * sqrt(2)/sqrt(pi);
                [sigma, ci] = mle(xx, 'distribution', 'Half Normal');
                
                fprintf('%s: halfnorm mean, [95%%CI] = %0.3f°, [%0.3f°, %0.3f°]\n', ...
                    str, rad2deg(sig2mean(sigma)), ...
                    rad2deg(sig2mean(ci)));
        end
        
    end

   
   function T = get_iw_table(min_electrodes)
       if nargin < 1, min_electrodes = 0; end
        T = readtable('iw_table');
        T(T.nchannels < min_electrodes, :) = [];
        
        T.seizure_dur = T.onset ./ T.onset_rel;
        T.phi_alt_deg = T.phi_alt/pi * 180;
        T.crossing_time = T.prop_sp;
        T.prop_sp = 10 * sqrt(2) * .4 ./ T.prop_sp;
        
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
    
   function [sc] = gscatter_pat(ax, xx, yy, pat, varargin)
        %[sc] = gscatter_pat(ax=gca, xx=[], yy, pat, ::'range', (Nx2)::)
        % only 6 because MG49 only shows up once now and this way colors
        % match with fig0
        cmap = repmat(lines(6), 2, 1);
        lso = 'oooooo^^^^^';
        
        switch nargin
            
            case 3 
                pat = [];
                if ~isa(ax, 'matlab.graphics.axis.Axes') %(xx, yy, pat)
                    pat = yy;
                    yy = xx;
                    xx = ax;
                    ax = gca;
                end
                if ~isnumeric(yy)  %(yy, pat, varargin) | (ax, yy, pat)
                    varargin = [pat varargin];
                    pat = yy;
                    yy = xx;
                    xx = [];
                end
            case 2 % (yy, pat) <- (ax, xx)
                pat = xx;
                yy = ax;
                xx = [];
                ax = gca;
            otherwise
                if ~isa(ax, 'matlab.graphics.axis.Axes')  % (xx, yy, pat, varargin)
                    varargin = [pat varargin];
                    pat = yy;
                    yy = xx;
                    xx = ax;
                    ax = gca;
                end
                
                if ~isnumeric(yy)  % (yy, pat, varargin) | (ax, yy, pat, varargin)
                    varargin = [pat varargin];
                    pat = yy;
                    yy = xx;
                    xx = [];
                end
                
        end
        assert(iscell(pat) || isstring(pat), ...
            'Expected F.(ax*, xx*, yy, pat). Asterisks denote optional.');
        sts = ax.NextPlot;
        
        range_ind = false(size(varargin));
        range_ind(1:2:end) = contains(lower(varargin(1:2:end)), 'range');
        if any(range_ind)
            data_range = varargin{find(range_ind) + 1};
        else
            data_range = nan(numel(yy), 2);
        end
        assert(size(data_range, 1) == numel(yy));
        
        
        % Create mapping and reorder data
        sz = SeizureInfo;
        [~, p, pa] = findgroups(sz.patient, sz.patientAlt);
        MAP = containers.Map(p, pa);
        G2 = cellfun(@(x) MAP(x), pat);
        
        [G2, so] = sort(G2);
        yy = yy(so);
        pat = pat(so);
        data_range = data_range(so, :);
        
        % maintain the same colors/lines every time
        

        if isempty(xx)
            xx = rescale(G2, .8, 1.2, ...
                'inputmin', 1, 'inputmax', 11);
            XLIM = [.6 1.4];
            XTICK = [];
            if ~all(isnan(data_range))
                xx = splitapply(@(x) {x(1)+.15*normalize(1:numel(x), 'center')'}, G2, G2);
                xx = cat(1, xx{:});
                XLIM = quantile(xx, [0 1]) + .1*range(xx) * [-1 1];
            end
        else
            xx = xx(so);
            XLIM = [];
            XTICK = 'auto';
        end

        for ii = 1

            ln = plot(ax, [xx(:), xx(:)]', data_range');
            hold(ax, 'on');
            set(ax, 'colororder', cmap(G2, :));
            sc = gscatter(ax, xx, yy, pat, ...
                cmap(unique(G2), :), lso(unique(G2)), [], 'on');
            for ss = sc', ss.MarkerFaceColor = ss.Color; end
            hold(ax, 'off');
            grid on
            box off
            xticks(XTICK);
            xlabel([]);
            if ~isempty(XLIM), xlim(XLIM); end
            
        end
        ax.NextPlot = sts;

    end

end
    
    
end


%% --- Local functions ---

function save_fig_(F, h)

outname = F.prefix(h.Tag);
F.print(h, outname); 
% savefig(h, outname)


end




