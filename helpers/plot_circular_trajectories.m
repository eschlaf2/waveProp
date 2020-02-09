function plot_circular_trajectories(res, title_str)
		
gcf; fullwidth(true);
if ~exist('name', 'var'), title_str = ''; end
nF = numel(res);
ax = gobjects(2*nF, 1);
for ii = 1:nF
	fields = fieldnames(res(ii).data);
	ax(ii) = subplot(2, nF, ii, polaraxes());
	ax(ii + nF) = subplot(2, nF, ii + nF, polaraxes());
	[ax(ii), ax(ii + nF)] = ...
		plot_wave_polar(res(ii), ax(ii), ax(ii + nF));

end
legend(rename_metrics(fields), 'position', [.9 .55 0 0])
rlim = arrayfun(@(a) a.RLim(2), ax(1:nF));
for kk = 1:numel(ax)
	if kk <= nF, ax(kk).RLim = [0 max(rlim)]; end
	ax(kk).ThetaTickLabel = [];
	ax(kk).RTickLabel = [];
end
ttl = @(s) annotation('textbox', ...
	'string', s, ...
	'Position', [.5 .05 0 0], ...
	'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
	'FitBoxToText', 'on', ...
	'LineStyle', 'none');
ttl(title_str)

end