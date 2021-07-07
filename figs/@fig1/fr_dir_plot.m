function fr_dir_plot(F, sub)

assert(ismember(sub, 'BCDEH'))
mm = F.Metrics.(sub);
mea = MEA(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', ...
    F.SimType, F.SimType, F.DetailSim));
style = F.Style;
fits = WaveProp.load(mea, {mm});
fits.(mm).MinFinite = 1;


ax = axes(figure);
plot(ax, mea.AllTime, rescale(mean(mea.firing_rate, 2), -pi, pi), 'color', .75*[1 1 1], 'linewidth', .5);
hold(ax, 'on');
scatter(ax, fits.(mm).time, fits.(mm).Direction, 'filled', 'markerfacecolor', style.(mm).color, 'sizedata', 20);
hold(ax, 'off');
axis(ax, 'tight');
ylim(ax, [-pi pi]);
yticks(ax, [-pi 0 pi]);
yticklabels(ax, {'-\pi', '0', '\pi'});
ax.Tag = mm;
h = ax.Parent;
h.Tag = mm;
F.(sub) = h;

end




