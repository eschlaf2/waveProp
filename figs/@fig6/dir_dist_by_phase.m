function dir_dist_by_phase(F, metrics, class)
    % dir_dist_by_phase(F, metrics=F.Metrics, class=1)
    % 
    
    if nargin < 3 || isempty(class), class = inf; end
    if nargin < 2 || isempty(metrics), metrics = F.Metrics; end
    if ischar(metrics), metrics = {metrics}; end
    
    
    % set some constants
    data_all = F.CompiledData;
    data_all(data_all.angle_p_val > .05, :) = [];
    data_all.dir = angle(exp(1j*(data_all.dir - data_all.iw_angle)));
    switch class
        case 1
            idx = strcmpi('u', data_all.ubo);
%                     idx = data_all.phase == 3;
        case 2
            idx = strcmpi('b', data_all.ubo);
%                     idx = data_all.phase == 1;
        case Inf
            idx = ismember(lower(data_all.ubo), {'u', 'b', 'o'});
            
        otherwise
            error('Groupname not recognized')
    end
    
    
    GROUPNAMES = ["Post" "Pre"];
    AXW = 2;  % width of single axis
    AXH = .8 * AXW;
    LMARGIN = 0.4;
    BMARGIN = 0.25;
    SPACING = 0.35;
    NMET = numel(metrics);
    FW = AXW * NMET + SPACING * (NMET-1) + 2 * LMARGIN;  % figwidth
    LS = {'-', ':'};
    
    
    
    % create fig and axes
    h = figure('units', 'inches', 'position', [0 0 FW AXH+LMARGIN]); % , ...
%         'resize', 'off');
    ax = gobjects(NMET, 1);
    for ii = 1:NMET
        left = LMARGIN + (AXW + SPACING) * (ii-1);
        ax(ii) = axes(h, 'units', 'inches', ...
            'position', [left BMARGIN AXW AXH]); 
    end
    
    for ii = 1:NMET
        mm = metrics{ii};
        color = [F.Style.(mm).color; .25 * [1 1 1]];
%         patch_colr = [.5 * [1 1 1]; color];
        
        % show each phase
        for jj = 2:-1:1
            switch GROUPNAMES(jj)
                case "Pre"
                    t_mask = data_all.phase == 1;
                case "Post"
                    t_mask = data_all.phase == 3;
                otherwise
                    error('Groupname not recognized')
            end
            data = data_all(idx & t_mask, :);
            

            % get distributions for all mm analyses
            d0 = data(strcmpi(data.metric, mm), :);
%             [G, pat, sz] = findgroups(d0.patient, d0.seizure); %#ok<ASGLU>
            [G, pat] = findgroups(d0.patient); %#ok<ASGLU>
            gN = splitapply(@numfinite, d0.dir, G);
            mask = ismember(G, find(gN < 50));
            d0(mask, :) = [];
            [G, pat] = findgroups(d0.patient); %#ok<ASGLU>
            % Bootstrap
            % Compute a distribution comprised of N samples from each 
            % patient. Repeat this BN times to generate stats (mean/ci) of
            % this distribution.
            xi = linspace(-pi/2, 3*pi/2, 128);
            [mn, ci] = bootstrap_(d0, G, xi);
            

            % compute histograms
            dens = ...
                splitapply(@(x) circ_ksdens(x, xi), d0.dir, G);  


            % compute mn and ci *** update this! ***
%             mn = mean(dens);
%             sd = std(dens);
%             ci = 1.96*sd ./ sqrt(size(dens, 1));
%             [mn, ci] = bootstrap_(dens, xi, N);


            % create patch
            xx = xi' / pi * 180;
%             low = [xx, mn - ci];
%             high = [xx, mn + ci];
            low = [xx, ci(:, 1)];
            high = [xx, ci(:, 2)];
            faces = 1:2*numel(xx);
            verts = [low; flipud(high)];
            patch(ax(ii), 'faces', faces, 'vertices', verts, ...
                'facecolor', color(jj, :), 'facealpha', .3, ...
                'edgecolor', 'none', 'tag', 'ci');

            hold(ax(ii), 'on');
            plot(ax(ii), xx, dens', 'linewidth', .75, ...
                'color', color(jj, :), 'tag', 'individual');

            xticks(ax(ii), -90:90:270)
            yticks(ax(ii), 0:.1:1.0);
            %     legend(pat(idx));

            plot(ax(ii), xx, mn, LS{jj}, 'color', color(jj, :), ...
                'linewidth', 2, 'displayname', GROUPNAMES(jj), ...
                'tag', 'mean')
        end
        hold(ax(ii), 'off');
        grid(ax(ii), 'on');
        title(ax(ii), F.MetricNames.(mm), 'color', F.Style.(mm).color);
        axis(ax(ii), 'tight');
        lgd = findobj(ax(ii), 'tag', 'mean');
        legend(ax(ii), lgd, 'location', 'southeast');


        % reorder children
        tags = arrayfun(@(a) a.Tag, ax(ii).Children, 'uni', 0);
        individs = strcmpi(tags, 'individual');
%         mns = find(strcmpi(tags, 'mean'));
%         ci = find(strcmpi(tags, 'ci'));
        set(ax(ii).Children(individs), 'visible', 'off');
%         ax(ii).Children = ax(ii).Children([mns; individs; ci]);
        
    end
    yy = get(ax, 'YLim');
    yy = quantile(cat(1, yy{:}), [0 1], 'all');
    set(ax, 'xlim', [-90 270]);
    set(ax, 'ylim', [0 yy(2)]);
%     set(ax, 'xlim', [-180 180]);
    F.print(h, F.prefix(['UvB_class' num2str(class)]));
end


%% local functions

function [mn, ci, dist, samples] = bootstrap_(d0, G, xi) 
% [mn, ci, dist, samples] = bootstrap_(d0, G, xi)
% Compute a distribution comprised of N samples from each 
% patient. Repeat this process BN times to generate stats (mean/ci) of
% this distribution.    
    
    
    % need d0 and G
    N = 500;
    BN = 1000;
    
    % compute std for each patient to add noise
    sd = splitapply(@(x) circ_std(x, [], [], 'omitnan'), d0.dir, G); 
%     sd = sd * .5; 
    samples = nan(N, max(G), BN);
    dist = nan(128, BN);
    
    for bb = 1:BN
        for ii = 1:max(G)
            dat = d0.dir(G == ii);
            inds = randi(numel(dat), N, 1);  % random indices
            samples(:, ii, bb) = dat(inds) + randn(size(inds)) * sd(ii);
            inds = randi(numel(dat), N, 1);  % random indices
            samples(:, ii, bb) = dat(inds) + randn(size(inds)) * sd(ii);
        end
        dist(:, bb) = circ_ksdens(reshape(samples(:, :, bb), [], 1), xi);
    end

    samples = reshape(samples, [], BN);
    ci = quantile(dist, [.025 .975], 2);
    mn = mean(dist, 2);


end


function [mn, ci] = bootstrap_FINAL(dens) 
    N = 1e3;
    Ns = size(dens, 1);
    mn = zeros(N, size(dens, 2));
    for ii = 1:N
        inds = randi(Ns, [Ns, 1]);
        mn(ii, :) = mean(dens(inds, :)); 
    end
    
    ci = quantile(mn, [.025 .975]);
    mn = mean(mn);
%     inds = randi(numel(dirs), N);
%     boot = dirs(inds);
end

function [dens, xi] = compute_hist_FINAL(dirs)
    xi = linspace(-pi/2, 3*pi/2, 128);
    d0 = dirs + [2*pi 0 -2*pi];
    dens = ksdensity(d0(:), xi, 'bandwidth', pi/16);
    dens = dens/sum(dens) * 2*pi;
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

