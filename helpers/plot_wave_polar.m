function [res, ax1, ax2] = plot_wave_polar(res, sig, ax1, ax2)
    
    if ~exist('ax1', 'var'); ax1 = polaraxes(); end
    if ~exist('ax2', 'var'); ax2 = polaraxes(); end
    
    fields = fieldnames(res.data);
% 	whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
%     res.time = res.data.(fields{whichfields(1)}).computeTimes / 1e3;
	time = res.time - min(res.time);
	[res.Z, res.Vx, res.Vy] = ...
		deal(zeros(length(time), length(fields)));
	for jj = 1:numel(fields)
		for f = 'Zp'
			
			switch f
				case {'Z', 'p'}
					data = interp1(...
						res.data.(fields{jj}).computeTimes / 1e3, ...
						res.data.(fields{jj}).(f), time, 'nearest');
				case 'Vx'
					data = interp1(...
						res.data.(fields{jj}).computeTimes / 1e3, ...
						res.data.(fields{jj}).V(1, :), time, 'nearest');
				case 'Vy'
					data = interp1(...
						res.data.(fields{jj}).computeTimes / 1e3, ...
						res.data.(fields{jj}).V(2, :), time, 'nearest');
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
%     Zu = unwrap(res.Z);
%     valid = ~isnan(Zu);
%     f = fit(res.time(valid)', Zu(valid), 'smoothingspline', 'smoothing', .05);
	data1 = smoothdata(...
		interp1(res.time, exp(1j*res.Z), tt), ...
		'movmean', 5 / (tt(2) - tt(1)), 'includenan');
    
	polarplot(ax1, angle(data1), tt, '-', 'linewidth', 2);
	axis tight;
    title(strrep(res.name, '_', ''));
    
% 	data2 = cumsum(interp1(res.time, exp(1j*res.Z), tt), 'omitnan');
	data2 = cumsum(exp(1j*res.Z), 'omitnan');
	polarplot(ax2, angle(data2), time, '.', ...
		'markersize', 10);
	hold off;
	axis tight
    title(strrep(res.name, '_', ''));
end