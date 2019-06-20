% pat = 'MG49';
files = dir([pat '_Seizure*wave_prop_all_waves.mat']);
nF = numel(files);

clear res;
res(nF) = struct(...
    'name', [], ...
    'data', [], ...
    'Z', [], ...
    'time', [], ...
    'Vx', [], ...
    'Vy', []);

figure(1); clf; fullwidth(true);
whichfields = [1:3 5:6]; %[2 3 4];
whichfields = 1;
for ii = 1:nF
	res(ii).name = files(ii).name(strfind(files(ii).name, 'Seizure')+(7:8));
	res(ii).data = load(fullfile(files(ii).folder, files(ii).name));
	fields = fieldnames(res(ii).data);
	res(ii).time = res(ii).data.(fields{whichfields(1)}).computeTimes / 1e3;
	[res(ii).Z, res(ii).Vx, res(ii).Vy] = ...
		deal(zeros(length(res(ii).data.(fields{whichfields(1)}).computeTimes), length(fields)));
	for jj = 1:numel(fields)
		for f = {'Z', 'Vx', 'Vy'}
			f = f{1};
			switch f
				case 'Z'
					data = res(ii).data.(fields{jj}).(f);
				case 'Vx'
					data = res(ii).data.(fields{jj}).V(1, :);
				case 'Vy'
					data = res(ii).data.(fields{jj}).V(2, :);
			end
			res(ii).(f)(:, jj) = data;
            % Remove values where fit is not significant
			res(ii).(f)(res(ii).data.(fields{jj}).p > .05, :) = nan;
            % Remove values where slope is zero in both directions
            res(ii).(f)(all(abs(res(ii).data.(fields{jj}).beta(1:2, :)) < eps), :) = nan;
		end
	end
	ax(ii) = polaraxes();
% 	ax(ii) = axes();
	subplot(2, nF, ii, ax(ii));
	tt = linspace(res(ii).time(1), res(ii).time(end), length(res(ii).time) * 10);
% 	data = interp1(res(ii).time, smoothdata(unwrap(res(ii).Z(:, whichvars)), 'movmean', 20), tt);
    Zu = unwrap(res(ii).Z(:, whichfields));
%     valid = ~isnan(Zu);
%     f = fit(res.time(valid)', Zu(valid), 'smoothingspline', 'smoothing', .05);
	data = smoothdata(...
		interp1(res(ii).time, Zu, tt), ...
		'rlowess', 1/mean(diff(res(ii).time))^3, 'omitnan');
	polarplot(ax(ii), data, tt, '-', 'linewidth', 2);
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
	polarplot(ax(ii + nF), res(ii).Z(:, whichfields), res(ii).time, '.', ...
		'markersize', 10);
	hold off;
	axis tight
    title(strrep(res(ii).name, '_', ''));
    
end
legend(fields(whichfields), 'position', [.9 .55 0 0])
rlim = arrayfun(@(a) a.RLim(2), ax);
for kk = 1:numel(ax)
	ax(kk).RLim = [0 max(rlim)];
	ax(kk).ThetaTickLabel = [];
	ax(kk).RTickLabel = [];
end
ttl = @(s) annotation('textbox', ...
    'string', s, ...
    'Position', [.5 .98 0 0], ...
    'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'FitBoxToText', 'on', ...
    'LineStyle', 'none');
ttl(pat)
%%
figure(2); clf; fullwidth(true);
S = 3;
metrics = [1 3:6];
pairs = nchoosek(metrics, 2);
nP = size(pairs, 1);
r = floor(sqrt(nP)); c = ceil(sqrt(nP));
fields = fieldnames(res(S).data);
for dd = 1:nP
    ax = polaraxes();
    subplot(r, c, dd, ax)
    d1 = pairs(dd, 1); d2 = pairs(dd, 2);
    data = res(S).Z(:, d1) - res(S).Z(:, d2);
    edges = linspace(0, 2*pi, 41);
    bc = histcounts(mod(data, 2*pi), edges);
    bc = bc/max(bc) * res(S).time(end);
    polarhistogram(ax, 'BinEdges', edges, 'BinCounts', bc); hold on;
    polarplot(ax, data, res(S).time, '.'); hold off;
    title(ax, sprintf('%s - %s', fields{d1}, fields{d2}))
    axis tight
    ax.ThetaTickLabel = [];
%     set(gca, 'thetaticklabels', [])
end
ttl(sprintf('%s Seizure%s', pat, strrep(res(S).name, '_', '')));
