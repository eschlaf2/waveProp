function [ax1, ax2] = plot_wave_polar(res, ax1, ax2)
    
    if ~exist('ax1', 'var'); ax1 = polaraxes(); end
    if ~exist('ax2', 'var'); ax2 = polaraxes(); end
   
% 	data = interp1(res.time, smoothdata(unwrap(res.Z(:, whichvars)), 'movmean', 20), tt);
%     Zu = unwrap(res.Z);
%     valid = ~isnan(Zu);
%     f = fit(res.time(valid)', Zu(valid), 'smoothingspline', 'smoothing', .05);
	tt = min(cat(1, res.time)):.01:max(cat(1, res.time));	
	data1 = smoothdata(...
		interp1(res.time, exp(1j*res.Z), tt), ...
		'movmean', .1 / (tt(2) - tt(1)), 'includenan');
    
	polarplot(ax1, angle(data1), tt - tt(1), '.', 'linewidth', 2, 'markersize', 4);
	axis tight;
    title(strrep(res.name, '_', ''));
    
	data2 = cumsum(interp1(res.time, diff(tt(1:2)) * exp(1j*res.Z), tt), 'omitnan');
% 	data2 = cumsum(exp(1j*res.Z), 'omitnan');
	polarplot(ax2, angle(data2), cumsum(tt), ...  % abs(data2), ...
                '.', 'linewidth', 2, ...
		'markersize', 4);
	hold off;
	axis tight
    title(strrep(res.name, '_', ''));
end
