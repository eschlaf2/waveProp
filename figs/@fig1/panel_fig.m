function h = panel_fig(F)
fname = sprintf('SCM/%s/%s_%d_info.mat', F.SimType, F.SimType, F.DetailSim);
load(fname, 'params', 'Qe_movie');
h = figure(); fullwidth;
set(h, 'units', 'points', 'position', [30   500   939   294])
set(h, 'defaultAxesFontSize', 13)

t = tiledlayout(h, 3, 10);
set(t, 'padding', 'compact', 'tilespacing', 'compact');

ctr = params.electrodes.centerNP;
ext = params.electrodes.dimsNP;
% stim = params.model.stim_center;
stim = nan(size(params.source, 3), 2);
for ii = 1:size(params.source, 3)
    [xi, xj] = find(params.source(:, :, ii) == max(params.source(:, :, ii), [], 'all'));
    stim(ii, :) = mean([xi xj]);
end
x = ctr(1) + [1 ext(1) ext(1) 1] - round(ext(1)/2);
y = ctr(2) + [1 1 ext(2) ext(2)] - round(ext(2)/2);
np =@(ax) patch(ax, x, y, 'red', 'linestyle', 'none', 'facealpha', .5);
for ii = 11:40
    ax = nexttile(t, ii-10); 
    contourf(ax, Qe_movie(ii).cdata, 'linestyle', 'none'); 
    colormap(ax, 1 - gray); 
    np(ax); 
    hold(ax, 'on')
    plot(ax, stim(:, 2), stim(:, 1), 'r.', 'markersize', 15);
    hold(ax, 'off')
    title(num2str(ii-10)); 
end
axis(t.Children, 'square');
set(t.Children, 'xtick', [], 'ytick', []);
% set(h, 'position', [0    0.5633    0.6521    0.3222])
h.Tag = 'scm_panel';
F.A = h;
end

