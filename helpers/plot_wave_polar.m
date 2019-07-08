function [ax1, ax2] = plot_wave_polar(res, metrics, sig)
    
    fields = fieldnames(res.data);
	whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
    res.time = res.data.(fields{whichfields(1)}).computeTimes / 1e3;
	time = res.time - min(res.time);
	[res.Z, res.Vx, res.Vy] = ...
		deal(zeros(length(res.data.(fields{whichfields(1)}).computeTimes), length(fields)));
	for jj = 1:numel(fields)
		for f = 'Zp'
			switch f
				case {'Z', 'p'}
					data = res.data.(fields{jj}).(f);
				case 'Vx'
					data = res.data.(fields{jj}).V(1, :);
				case 'Vy'
					data = res.data.(fields{jj}).V(2, :);
			end
			res.(f)(:, jj) = data;
            % Remove values where fit is not significant
			mask = res.data.(fields{jj}).p < sig;
			res.(f)(~mask, jj) = nan;
            % Remove values where slope is zero in both directions
            res.(f)(all(abs(res.data.(fields{jj}).beta(1:2, :)) < eps), jj) = nan;
		end
	end
% 	ax = axes();
	tt = (time(1): 1e-3: time(end));
% 	data = interp1(res.time, smoothdata(unwrap(res.Z(:, whichvars)), 'movmean', 20), tt);
    Zu = unwrap(res.Z(:, whichfields));
%     valid = ~isnan(Zu);
%     f = fit(res.time(valid)', Zu(valid), 'smoothingspline', 'smoothing', .05);
	data = smoothdata(...
		interp1(res.time, Zu, tt), ...
		'movmean', 5 / (tt(2) - tt(1)), 'omitnan');
    
    ax1 = polaraxes();
	polarplot(ax1, data, tt, '-', 'linewidth', 2);
	axis tight;
    title(strrep(res.name, '_', ''));
    
    ax2 = polaraxes();
	polarplot(ax2, res.Z(:, whichfields), time, '.', ...
		'markersize', 10);
	hold off;
	axis tight
    title(strrep(res.name, '_', ''));
end