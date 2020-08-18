function [H, centers, T] = hist(obj, window, mean_center, show_directions, h)
    % [H, edges, T] = hist(window=1, mean_center=true, show_directions=true, h=gcf)

    % Parse inputs
    if nargin < 5, h = gcf; end
    if nargin < 4 || isempty(show_directions), show_directions = true; end
    if nargin < 3 || isempty(mean_center), mean_center = true; end
    if nargin < 2 || isempty(window), window = 1; end
    if mean_center
        obj.Direction = angle(exp(1j*(obj.Direction - obj.mean)));
    end
    TL = tiledlayout(h, 1, 9);
    ax_summ = nexttile(TL, [1, 2]);
    ax_main = nexttile(TL, [1, 7]);
    EDGES = linspace(-pi, pi, 65);
    
    % Compute windowed empirical distributions
    [H, T] = compute_moving_hist_(obj, EDGES, window);
    
    % Smooth result 
    [H, centers] = smooth_circular_(H, T, EDGES, window);
    
    % Create 2D contour plot of windowed histograms
    hist2D_(ax_main, T, centers, H, show_directions);

    % Create 1D plot of direction histograms
    hist1D_(ax_summ, obj, EDGES);

    % Tag axes so they are easy to identify later
    ax_main.Tag = 'Main';
    ax_summ.Tag = 'Summ';

end


%% Local functions
function [H, T] = compute_moving_hist_(obj, edges, window)
    t = obj.time;
    T = t(1):.01:t(end)+.01;
    H = nan(length(edges)-1, length(T));
    for ii = 1:length(T)
        t_inds = abs(t - T(ii)) <= window/2;
        dir = obj.Direction(t_inds);
        H(:, ii) = histcounts(dir, edges) ./ window;
    end
end

function [H, centers] = smooth_circular_(H, T, edges, window)
    % Performs smoothing of directional data by replicating the data in the
    % +/- 2pi directions, smoothing and then selecting the data in interval
    % [-pi pi]
    win = gausswin(round(pi/4 / diff(edges(1:2)))) ...
        .* gausswin(round(window ./ diff(T(1:2))))';
    win = win / sum(win, 'all');
    H = [H; H; H];
    centers = edges(1:end-1)' + diff(edges(1:2))/2;
    centers = [centers - 2*pi; centers; centers + 2*pi];
    H = conv2(H, win, 'same');
    angle_mask = abs(centers) <= 1.25*pi;
    centers = centers(angle_mask);
    H = H(angle_mask, :);
end

function hist2D_(ax, T, centers, H, show_directions)
    contourf(ax, T, centers, H, quantile(H(:), .1:.05:.9), ...
        'linestyle', 'none'); 
    if show_directions
        NP = ax.NextPlot;
        ax.NextPlot = 'add';
        plot(ax, obj.time, obj.Direction, '.', ...
            'markersize', 10);
        ax.NextPlot = NP;
    end
    colorbar(ax);
    xlabel(ax, 'Time (s)')
    ylabel(ax, 'Direction')
    yticks(ax, [-pi, pi]);
    grid(ax, 'on');
end

function hist1D_(ax, obj, edges)
    H = histcounts(obj.Direction, edges, 'normalization', 'pdf');
    centers = movmean(edges, 2); centers = centers(2:end);
    H = smoothdata([H, H, H]);
    C = [centers - 2*pi, centers, centers + 2*pi];
    mask = abs(C) <= 1.25*pi;
    plot(ax, H(mask), C(mask), ...
        'linewidth', 2);
    axis(ax, 'tight');
    xticks(ax, []); yticks(ax, [-pi pi]); grid(ax, 'on');
    
end

