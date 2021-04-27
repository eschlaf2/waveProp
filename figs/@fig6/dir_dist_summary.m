function dir_dist_summary(F, metrics, period)
    
    if nargin < 3 || isempty(period), period = ''; end
    if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
    if ischar(metrics), metrics = {metrics}; end
    
    
    % set some constants
    data_all = F.CompiledData;
    data_all.dir = angle(exp(1j*(data_all.dir - data_all.iw_angle)));
    
    switch lower(period)
        case 'early', t_mask = data_all.phase == 1;
        case 'late', t_mask = data_all.phase == 3;
        otherwise, t_mask = true(size(data_all.phase));
    end
    GROUPNAMES = {'Class I', 'Class II'};
    AXW = 1.5;  % width of single axis
    AXH = .8 * AXW;
    LMARGIN = 0.4;
    BMARGIN = 0.25;
    SPACING = 0.35;
    FW = AXW * 3 + SPACING * 2 + 2 * LMARGIN;  % figwidth
    NMET = numel(metrics);
    LS = {'-', ':'};
    
    
    
    % create fig and axes
    h = figure('units', 'inches', 'position', [0 0 FW AXH+LMARGIN], 'resize', 'off');
    ax = gobjects(3, 1);
    for ii = 1:NMET
        left = LMARGIN + (AXW + SPACING) * (ii-1);
        ax(ii) = axes(h, 'units', 'inches', ...
            'position', [left BMARGIN AXW AXH]); 
    end
    
    for ii = 1:NMET
        mm = metrics{ii};
        color = [F.Style.(mm).color; .25 * [1 1 1]];
%         patch_colr = [.5 * [1 1 1]; color];
        
        % show each class
        for jj = 2:-1:1
            switch GROUPNAMES{jj}
                case 'Class I'
                    idx = strcmpi('u', data_all.ubo);
%                     idx = data_all.phase == 3;
                case 'Class II'
                    idx = strcmpi('b', data_all.ubo);
%                     idx = data_all.phase == 1;
                otherwise
                    error('Groupname not recognized')
            end
            data = data_all(idx & t_mask, :);
            

            % get distributions for all mm analyses
            d0 = data(strcmpi(data.metric, mm), :);
            [G, pat, sz] = findgroups(d0.patient, d0.seizure); %#ok<ASGLU>

            % compute histograms
%             [hist, centers, pks, locs] = ...
            [dens, xi] = ...
                splitapply(@compute_hist_, d0.dir, G);  % ** local fun **


            % compute mn and ci *** update this! ***
            mn = mean(dens);
            sd = std(dens);
            ci = 1.96*sd ./ sqrt(size(dens, 1));


            % create patch
            xx = xi(1, :) / pi * 180;
            low = [xx; mn - ci]';
            high = [xx; mn + ci]';
            faces = 1:2*numel(xx);
            verts = [low; flipud(high)];
            patch(ax(ii), 'faces', faces, 'vertices', verts, ...
                'facecolor', color(jj, :), 'facealpha', .3, ...
                'edgecolor', 'none', 'tag', 'ci');

            hold(ax(ii), 'on');
            plot(ax(ii), xx, dens', 'linewidth', .75, ...
                'color', color(jj, :), 'tag', 'individual');

            xticks(ax(ii), -90:90:270)
            yticks(ax(ii), 0:.1:1);
            %     legend(pat(idx));

            plot(ax(ii), xx, mn, LS{jj}, 'color', color(jj, :), ...
                'linewidth', 2, 'displayname', GROUPNAMES{jj}, ...
                'tag', 'mean')
        end
        hold(ax(ii), 'off');
        grid(ax(ii), 'on');
        title(ax(ii), F.MetricNames.(mm), 'color', F.Style.(mm).color);
        axis(ax(ii), 'tight');
        lgd = findobj(ax(ii), 'tag', 'mean');
        legend(ax(ii), lgd);


        % reorder children
        tags = arrayfun(@(a) a.Tag, ax(ii).Children, 'uni', 0);
        individs = strcmpi(tags, 'individual');
%         mns = find(strcmpi(tags, 'mean'));
%         ci = find(strcmpi(tags, 'ci'));
        set(ax(ii).Children(individs), 'visible', 'off');
%         ax(ii).Children = ax(ii).Children([mns; individs; ci]);
        
    end
%     linkaxes(ax, 'y');
    set(ax, 'xlim', [-90 270]);
%     set(ax, 'xlim', [-180 180]);
    F.print(h, F.prefix(['UvB' period]));
end


%% local functions

function boot = bootstrap_(hist) %#ok<INUSD,DEFNU>
    N = 1e4;
    inds = randi(numel(dirs), N);
    boot = dirs(inds);
end

function [dens, xi] = compute_hist_(dirs)
    xi = linspace(-pi/2, 3*pi/2, 128);
    d0 = dirs + [2*pi 0 -2*pi];
    dens = ksdensity(d0(:), xi, 'bandwidth', pi/16);
end

function [hist, centers, pks, locs] = ZZcompute_hist_(dirs)
    CENTERS = linspace(-pi, pi, 129);
    binwidth = diff(CENTERS(1:2));
    edges = CENTERS - binwidth/2;
    counts = histcounts(dirs, edges, 'normalization', 'pdf');
    
    centers = CENTERS(1:end-1);
    
    cc = [centers-2*pi centers centers + 2*pi];
    counts_temp = [counts counts counts];
    counts_sm = smoothdata(counts_temp, 'gaussian', pi/4, 'SamplePoints', cc);
    
    % mask to original interval
    mask = abs(cc) <= pi;
    
    % normalize
%     counts_sm = counts_sm / max(counts_sm(mask));
    
    % find peaks
    [pks, locs] = findpeaks(counts_sm(mask), cc(mask), ...
        'NPeaks', 2, 'SortStr', 'descend', 'minpeakdistance', pi/4);
    
    cc = cc - locs(1) - pi/2;
    
    hist = interp1(cc, counts_sm, CENTERS, 'nearest');
    centers = (CENTERS / pi * 180) + 90;
        
end

