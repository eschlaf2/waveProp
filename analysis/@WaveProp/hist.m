function [H, centers, T] = ...
        hist(obj, window, mean_center, show_directions, h, angular_res, use_particle_smoother)
    % [H, edges, T] =  ...
    %       hist(window=0, mean_center=false, show_directions=true, ...
    %           h=gcf, angular_res=64, use_particle_smoother=false)

    % Parse inputs
    if nargin < 7, use_particle_smoother = false; end
    if nargin < 6 || isempty(angular_res), angular_res = 64; end
    if nargin < 5 || isempty(h), h = gcf; end
    if nargin < 4 || isempty(show_directions), show_directions = false; end
    if nargin < 3 || isempty(mean_center), mean_center = false; end
    if nargin < 2 || isempty(window), window = 1; end
    if mean_center
        obj.Direction = angle(exp(1j*(obj.Direction - obj.mean)));
    end
    if all(isnan(obj.Direction)), disp('No detections.'); [H, centers, T] = deal(nan); return; end
    TL = tiledlayout(h, 1, 9);
    ax_summ = nexttile(TL, [1, 2]);
    ax_main = nexttile(TL, [1, 7]);
%     EDGES = linspace(-pi, pi, angular_res + 1);
    xi = 1.25 * pi * linspace(-1, 1, round(1.5 * angular_res));
    
    % Compute windowed empirical distributions
    if ~use_particle_smoother
        [H, T, centers] = compute_moving_histn_(obj, window, xi);
    end
    
    
    % Create 2D contour plot of windowed histograms
    if use_particle_smoother
        obj.particle_smoother(ax_main, [], [], [], [], window);
    else
        hist2D_(ax_main, T, centers, H);
        
        if show_directions
            NP = ax_main.NextPlot;
            ax.NextPlot = 'add';
            plot(ax, obj.time, obj.Direction, '.', ...
                'markersize', 10);
            ax_main.NextPlot = NP;
        end
    end

    % Create 1D plot of direction histograms
    hist1D_(ax_summ, obj);

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

function [dens, sub_t, xi] = compute_moving_histn_(obj, window, xi)
    t = obj.time;
    sub_t = t(1):.1:t(end);
    dens = zeros(numel(xi), numel(sub_t));
    for ii = 1:numel(sub_t)
        inds = t > sub_t(ii) - window/2 & t < sub_t(ii) + window/2;
        if obj.HalfWin < .03, 
            dirs = obj.discharge_directions;
            dirs = dirs(inds);
        else
            dirs = obj.Direction(inds);
        end
        if sum(isfinite(dirs)) < 1, continue; end
        dirs = [dirs + 2*pi; dirs; dirs - 2*pi];
        dens(:, ii) = ksdensity(dirs, xi, 'bandwidth', pi/16); 
    end
end

function hist2D_(ax, T, centers, H)
    contourf(ax, T, centers, H, quantile(H(:), .1:.05:.9), ...
        'linestyle', 'none'); 
    
    colorbar(ax);
    xlabel(ax, 'Time (s)')
    ylabel(ax, 'Direction')
    yticks(ax, [-pi, pi]);
    grid(ax, 'on');
end

function hist1D_(ax, obj)
    dir = obj.Direction;
    [H, C] = ksdensity([dir - 2*pi; dir; dir + 2*pi], ...
        1.25 * pi * linspace(-1, 1, 100), ...
        'bandwidth', pi/16);
    plot(ax, H, C, 'linewidth', 2);
    axis(ax, 'tight');
    xticks(ax, []); yticks(ax, [-pi pi]); grid(ax, 'on');
    
end

