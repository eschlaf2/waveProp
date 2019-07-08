function plot_wave_polar(res, metrics, sig)
    
    fields = fieldnames(res(ii).data);
	whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
    res(ii).time = res(ii).data.(fields{whichfields(1)}).computeTimes / 1e3;
	time = res(ii).time - min(res(ii).time);
	[res(ii).Z, res(ii).Vx, res(ii).Vy] = ...
		deal(zeros(length(res(ii).data.(fields{whichfields(1)}).computeTimes), length(fields)));
	for jj = 1:numel(fields)
		for f = 'Zp'
			switch f
				case {'Z', 'p'}
					data = res(ii).data.(fields{jj}).(f);
				case 'Vx'
					data = res(ii).data.(fields{jj}).V(1, :);
				case 'Vy'
					data = res(ii).data.(fields{jj}).V(2, :);
			end
			res(ii).(f)(:, jj) = data;
            % Remove values where fit is not significant
			mask = res(ii).data.(fields{jj}).p < sig;
			res(ii).(f)(~mask, jj) = nan;
            % Remove values where slope is zero in both directions
            res(ii).(f)(all(abs(res(ii).data.(fields{jj}).beta(1:2, :)) < eps), jj) = nan;
		end
	end
	ax(ii) = polaraxes();
% 	ax(ii) = axes();
	subplot(2, nF, ii, ax(ii));
	tt = (time(1): 1e-3: time(end));
% 	data = interp1(res(ii).time, smoothdata(unwrap(res(ii).Z(:, whichvars)), 'movmean', 20), tt);
    Zu = unwrap(res(ii).Z(:, whichfields));
%     valid = ~isnan(Zu);
%     f = fit(res.time(valid)', Zu(valid), 'smoothingspline', 'smoothing', .05);
	data = smoothdata(...
		interp1(res(ii).time, Zu, tt), ...
		'movmean', 5 / (tt(2) - tt(1)), 'omitnan');
	polarplot(ax(ii), data, tt, '-', 'linewidth', 2);
    ax(ii).ColorOrder = varycolor(5);
% % 	plot3(ax(ii), cos(data), sin(data), tt); hold on
% 	ax(ii).ColorOrderIndex = 1;
% 	plot3(ax(ii), cos(res(ii).Z(:, whichfields)), ...
% 		sin(res(ii).Z(:, whichfields)), res(ii).time, '.', 'markersize', 10); hold off
% 	title([pat ' ' res(ii).name]);
	axis tight;
    title(strrep(res(ii).name, '_', '')); 
% 	hold on;
% 	set(gca, 'ColorOrderIndex', 1)
	ax(ii + nF) = polaraxes();
	subplot(2, nF, ii + nF, ax(ii + nF));
	polarplot(ax(ii + nF), res(ii).Z(:, whichfields), time, '.', ...
		'markersize', 10);
	hold off;
	axis tight
    title(strrep(res(ii).name, '_', ''));
end