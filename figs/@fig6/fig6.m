classdef fig6 < BVNY & handle

    
properties
%     Metrics = {'M', 'E', 'D1xwh'}
    Metrics = {'M', 'D10'}
    Seizures
    H
    PatientByRow = false  % if true, display each row as a patient; else each column is a patient
    Smoothing = 1  % seconds
    CompiledData
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

            if 1
            h = figure('name', [M.Name '_clust']); fullwidth(true);
            M.cluster(h);
            F.print(h, F.prefix_better(h.Name));
            close(h);
            end
        end
    end
    function distance_scatters(F, sz_num)
        % Shows the distance (correlation distance) between discharge
        % patterns and IW patterns
        PTHRESH = 5e-3;
        sz = F.SeizureInfo;
        if nargin < 2, sz_num = 1:numel(sz.patient); end
        if isa(sz_num, 'MEA'), mea = sz_num; sz_num = 1; end
%         h = figure('name', 'distance_scatters'); fullwidth(true);
%         T = tiledlayout(h, 'flow');
        for ii = 1:numel(sz_num)
            pat = sz.patient{sz_num(ii)};
            seizure = sz.seizure(sz_num(ii));
            if ~exist('mea', 'var')
            mea = MEA(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', pat, pat, seizure));
            end
            fprintf('%d: Starting %s\n', sz_num(ii), mea.Name)
            out = mea.get_IW_templates([], [], 8);
            assignin('base', 'out', out);
            assignin('base', 'mea', mea);
            
            
            M = WaveProp.load(mea, {'M10'});
            M = M.M10;

            
            if isempty(out.template), continue; end  % if there are no IW, skip the rest
            
            try  % tired of re-running; catch error and rerun later
            h2 = figure('name', ['IW_' mea.Name]); fullwidth(true);
            T2 = tiledlayout(h2, 5, 4);
            
            nexttile(T2)
            [~, D] = M.distance_scatter(out, PTHRESH);
            dist_rho = D.dist_rho;
            dist_p = D.dist_p;
            dist_t = D.time;

            mask = isfinite(dist_rho) & dist_p(:, 1) < PTHRESH;
            title(strrep(M.Name, '_', ' '));
            
            nexttile(T2)
            histogram(dist_p(mask, 2));
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

            
            if sum(mask) > 0
                gscatter(dist_t(mask) - out.time(dist_p(mask, 2)), ...
                    -(1 - dist_rho(mask))/2 + dist_p(mask, 2), ...   % don't forget the negative here since you put this in axis ij
                    dist_p(mask, 2), [], [], 6, 'off')

                grid on
                hold on

                [G, gp] = findgroups(dist_p(mask, 2));
                hc = histcounts(categorical(G));
                if any(hc == 1)
                    mask = mask & ~ismember(dist_p(:, 2), find(hc == 1));
                    [G, gp] = findgroups(dist_p(mask, 2));
                end
                mvmn = splitapply(@(x, t) ...
                    interp1(t, movmean(x, 10, 'SamplePoints', t), M.time', 'linear'), ...
                    1-dist_rho(mask), dist_t(mask), G);
                mvmn = fillmissing(mvmn, 'constant', 0);
                plot(M.time - out.time(gp)', -mvmn'/2 + gp', 'r', 'linewidth', 1)
                hold off
                xline(0)
                ylim([0 9])
                axis ij
            end
            title(sprintf('\\rho v. time (p<5e%d)', log10(PTHRESH/5)))
            
                        
            nexttile(T2)
            violinplot(1-dist_rho(mask), dist_p(mask, 2));
            bp = findobj(gca, 'type', 'scatter');
            set(bp, 'MarkerFaceAlpha', 1, 'SizeData', 6);
            ylim([-1 1])
            title('\rho by template')
            
            nexttile(T2)
            violinplot(log10(dist_p(:, 1)/5), dist_p(:, 2));
            bp = findobj(gca, 'type', 'scatter');
            set(bp, 'MarkerFaceAlpha', 1, 'SizeData', 6);
            title('p=5eYY by template')
            
            nexttile(T2)
            vp = violinplot(dist_t(mask), dist_p(mask, 2));
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
            plot(0:100, quantile(log10(dist_p(:, 1)/5), linspace(0, 1, 101)));
            yline(log10(PTHRESH/5));
            grid on
            title('Percentiles log10(p/5)')
            
            
            nexttile(T2)  % Leave a blank so the templates are all on the same row
            plot(0:100, quantile(abs(1-dist_rho(dist_p(:, 1) < PTHRESH)), linspace(0, 1, 101)))
            grid on
            ylim([0 1])
            title(sprintf('Percentiles abs(\\rho), p < 5e%d', log10(PTHRESH/5)))
            
            
            nexttile(T2)
            ksdensity(1-dist_rho, linspace(-1, 1, 1000), 'bandwidth', .05);
            hold on; 
            ksdensity(1-dist_rho(mask), linspace(-1, 1, 1000), 'bandwidth', .05); 
            hold off;
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
                title({['Template ' num2str(jj)]; sprintf('t=%0.2fs', nanmedian(tpl, 'all'))})
                ylabel(sprintf('fr in [%d %d]', round(quantile(fr(inds), [0 1]))));
                xlim([0 11]);
                ylim([0 11]);
                xticklabels([]);
                yticklabels([]);
                colorbar;
                axis square
            end
            
            if isfinite(str2double(getenv('SGE_TASK_ID')))  % no print during testing
                F.print(h2, F.prefix_better(h2.Name));
                close(h2);
            end
            
            catch ME
                warning('error in %d: %s', sz_num(ii), mea.Name)
                disp(ME)
                for line = ME.stack'
                    fprintf('In %s (line %d)\n', line.name, line.line)
                end
%                 rethrow(ME)
                F.print(h2, F.prefix_better(h2.Name));
                close(h2);
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
            iw_info = load('iw_info.mat', fname);
            iw_info = iw_info.(fname);
            
            mea = MEA(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', pat, pat, seiz));
            iw_info_new = mea.get_IW_templates;
            
            for ff = F.Metrics
                W = fit.(ff{:});
                if strcmp(ff{:}, 'M10'), W.MinFinite = 10; end
                h = W.summary_plot;
                ax = findobj(h, 'Tag', 'direction_raster');
                for aa = findobj(h, 'type', 'axes')'
                    for t2 = iw_info_new.time(:)'
                        xline(aa, t2, 'r-', 'linewidth', 2);
                    end
                end
                if ~isempty(iw_info)
                    for aa = ax'
                        yline(aa, iw_info(4)/pi*180, 'k-', 'linewidth', 2);
                    end
                
                    for aa = findobj(h, 'type', 'axes')'
                        xline(aa, iw_info(1), 'k-', 'linewidth', 2);
                        
                    end
                
                pax = findobj(h, 'type', 'polar');
                hold(pax, 'on')
                polarplot(pax, [1; 1] * iw_info(4), [0; 1], 'k')
                hold(pax, 'off');
                end
                
                F.print(h, F.prefix_better(sprintf('%s_%s', fname, ff{:})));
%                 close(h)
            end
        end
    end
    
    function [out] = transient_iw(F, interval)
        % This will get all direction estimates within interval number of seconds before
        % and after the IW for each seizure and in each case return the mean direction
        % pre & post.
        if nargin < 2 || isempty(interval), interval = 5; end
        
        fits = F.compile_fits(interval);
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
        data = F.CompiledData;
        MTC = unique(data.metric);
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
            assert(all(strcmp(mtc, MTC)));  % throw an error if something is missing
            
            
            % Compute KL-divergence based on distribution of differences.
            % Compare dist of diffs to VM normal with same center and std
            W = WaveProp.load(sprintf('%s_Seizure%d', pat{ii}, sz(ii)), F.Metrics);
            for jj = 1:nP  
                fit1 = W.(mpairs{jj, 1});
                fit2 = W.(mpairs{jj, 2});
                if numel(fit1.time) > numel(fit2.time), 
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
        
        out = table(Dkl, Dmode, Dmean_diffs, Dstd_diffs, Dskew_diffs, Dcorr_pval, Dcorr_coeff, pat, sz);
        
    end
    function data = compile_fits(F, interval)
        % int is the interval surrounding the IW to assign phase 1 (early)
        % or phase 3 (late)
        
        if nargin < 2, interval = Inf; end  
        USE_MID = true;  % compute the intervals from the midpoint of the IW (or from the bounds)
        
        
        W = WaveProp.load([], F.Metrics);
        sz = F.SeizureInfo;
        data = table;

        for mm = string(F.Metrics)  % loop through metrics
            for ii = 1:height(sz)
                name = sprintf('%s_Seizure%d', sz.patient{ii}, sz.seizure(ii));
                iiW = cellfun(@any, strfind(W.Name, name));
                iw_info = load('iw_info.mat', name); iw_info = iw_info.(name);
                if isempty(iw_info), continue; end
    %             if isnan(iw_info(4)), iw_info(4) = 0; end
                fit = W.(mm)(iiW);
                fit.MinFinite = BVNY.MinFinite;
    %             fit.RotateBy = iw_info(3);
                data_temp = table(fit.Direction, fit.time, fit.AltMagnitude, 'VariableNames', {'dir', 'time', 'speed'});
                N = height(data_temp);
    %             name = strsplit(W.Name{ii}, {'_', 'Seizure'});
    %             data_temp.patient = repmat(string(name{1}), N, 1);
    %             data_temp.seizure = repmat(string(name{2}), N, 1);
                data_temp.patient = repmat(string(sz.patient{ii}), N, 1);
                data_temp.seizure = repmat(sz.seizure(ii), N, 1);
                data_temp.metric = repmat(mm, N, 1);
                data_temp.ubo = repmat(sz.UBO(ii), N, 1);
                
                phase = zeros(N, 1);
                iw0 = iw_info(2); % iw start
                iwF = iw_info(3); % iw end
                iwM = iw_info(1); % iw midpoint
                
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
                data_temp.iw_angle = repmat(iw_info(4), N, 1);
                data_temp.angle_p_val = repmat(iw_info(6), N, 1);
                data_temp.weight = 1/N * ones(N, 1);


                data = [data; data_temp];

            end
        end

    end

    function f = get_fname(F, ii)
        f = sprintf('%s_Seizure%d_fits.mat', ... 
            F.Seizures.patient{ii}, F.Seizures.seizure(ii));
    end
    function f = get_file(F, ii)
%        f = dir(F.get_fname(ii));
       f = F.get_fname(ii);
       f = f(1:end-9);
%        f = [F.Seizures.patient{ii} '_Seizure' num2str(F.Seizures.seizure(ii))];
    end
    
    function outname = prefix_better(F, tag)
        caller = dbstack(1);
        assignin('base', 'caller', caller);
        assert(isstring(tag) || ischar(tag));
        outname = sprintf('%sfig6/%s_%s_%d', F.Prefix, ...
            strrep(caller.name, '.', '_'), ...
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
        if nargin < 2 || isempty(patients), patients = unique(F.Seizures.patient); end
        if ischar(patients), patients = {patients}; end
        for pat = patients'
            try F.hist_figs(pat{:}, style); 
            catch ME, warning('make failed in %s', pat{:}); disp(ME); disp(ME.stack);
            end
        end
        F.dir_dist_by_phase();
%         F.dir_dist_by_phase([], 1);
%         F.dir_dist_by_phase([], 2);
%         F.dir_dist_summary([], 'early');
%         F.dir_dist_summary([], 'late');
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
    
    hist_figs(F, patient, style)
    hist_figs_speed(F, patient)
    dir_dist_summary(F, metrics, period)
    dir_dist_by_phase(F, metrics, class)
    plot_ttests(F)
end

methods (Static)
   [files, names] = txt2files(fname) 
   style = set_style
end
    
    
end


%% --- Local functions ---
function save_fig_(F, h)

outname = F.prefix(h.Tag);
F.print(h, outname); 
% savefig(h, outname)


end




