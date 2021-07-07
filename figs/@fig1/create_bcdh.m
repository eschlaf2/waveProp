function create_bcde(F, sub)

assert(ismember(sub, 'BCDH'))
mm = F.Metrics.(sub);
mea = MEA(sprintf('%s/%s_Seizure%d_Neuroport_10_10.mat', ...
    F.SimType, F.SimType, F.DetailSim));
style = F.Style;
fits = WaveProp.load(mea, 'metrics', mm);


ax = axes(figure);
plot(ax, mea.Time, rescale(mean(mea.firing_rate, 2), -pi, pi), 'color', .5*[1 1 1], 'linewidth', 1);
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




