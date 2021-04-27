function plot_ttests(F)
    % fig6 ttests
    h = figure('position', [0 0 3 1.6], 'resize', 'off');
    T = tiledlayout(h, 1, 3);
    ax = gobjects(numel(F.Metrics), 1);
    
    fits = WaveProp.load('metrics', F.Metrics);
    
    for mm = F.Metrics
        mag.(mm{:}) = arrayfun(@(f) nanmedian(f.AltMagnitude), fits.(mm{:}));
    end
    
    G = findgroups(F.SeizureInfo.UBO);

    G(G == 2) = nan;  % G = 2 => "other" (not uni or bi)
    G(G == 3) = 2;
    
    for ii = 1:numel(F.Metrics)
        ax(ii) = nexttile(T, ii);
        mm = F.Metrics{ii};
        [p, stats1, stats2] = boxplot_(mag.(mm), mm, ax(ii), G);
        display(mm, 'metric')
        display(p)
        stats1 %#ok<NOPRT>
        stats2 %#ok<NOPRT>
    end


    % print figure
    F.print(h, [F.Prefix 'fig6/fig6_ttests']);


end

%% local fun

function [p, stats1, stats2] = boxplot_(data, metric, aa, group)
    
%     BOX = [241,163,64]/255;
    BOX = [196,196,196; ...
        100, 100, 100]/255;
    MEDIAN = [1 1 1];
    OTHER = .2*[1 1 1];
    SCATTER = .2*[1 1 1];
    
%     disagree = sz.disagree > 1;
    mask = ~isnan(group);
    boxplot(aa, data(mask), group(mask));

    title(aa, metric);
    xticklabels(aa, {'B', 'U'});
    [~, p] = ttest2(data(group==1), data(group == 2));
    stats1 = fitdist(data(group==1), 'norm');
    stats2 = fitdist(data(group==2), 'norm');
    exponent = -log(p)/log(10);
    str = sprintf('p < 5e-%d', floor(exponent));

    hold(aa, 'on');
    scatter(aa, group, data, 25, SCATTER, 'filled')
    hold(aa, 'off');
    Y = aa.YLim;
    aa.YLim(2) = aa.YLim(2) + 0.25 * diff(Y);
    text(aa, 'string', str, ...
        'units', 'norm', 'position', [.95 .97 0], ...
        'horizontalalignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'fontweight', 'bold', ...
        'fontsize', 8);
    box_group = findobj(aa, 'tag', 'boxplot');
    set(box_group.Children, 'color', OTHER, 'linewidth', 2);
    box = findobj(aa, 'tag', 'Box');
    for ii = 1:numel(box)
        ln = box(ii);
        patch(ln.Parent, ln.XData, ln.YData, BOX(ii, :), ...
            'tag', 'Box', 'linestyle', 'none'); 
        delete(ln);
    end 
    box_group.Children = circshift(box_group.Children, -2);
    md = findobj(aa, 'tag', 'Median');
    set(md, 'color', MEDIAN);
end



