pat = 'MG49';
files = dir(['~/Desktop/temp/' pat '_Seizure*wave_prop.mat']);
nF = numel(files);

clear res;
res(nF) = struct('name', [], 'data', [], 'Z', [], 'time', [], 'Vx', 'Vy');
figure(); fullwidth(true);
whichvars = [2 3 4];
for ii = 1:nF
	res(ii).name = files(ii).name(strfind(files(ii).name, 'Seizure')+(7:8));
	res(ii).data = load(fullfile(files(ii).folder, files(ii).name));
	fields = fieldnames(res(ii).data);
	res(ii).time = res(ii).data.(fields{1}).computeTimes / 1e3;
	[res(ii).Z, res(ii).Vx, res(ii).Vy] = ...
		deal(zeros(length(res(ii).data.(fields{1}).computeTimes), length(fields)));
	for jj = 1:length(fields)
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
			res(ii).(f)(res(ii).data.(fields{jj}).p > .05, :) = nan;
		end
	end
% 	ax(ii) = polaraxes();
	ax(ii) = axes();
	subplot(2, nF, ii, ax(ii));
	tt = linspace(res(ii).time(1), res(ii).time(end), length(res(ii).time) * 10);
% 	data = interp1(res(ii).time, smoothdata(unwrap(res(ii).Z(:, whichvars)), 'movmean', 20), tt);
	data = smoothdata(...
		interp1(res(ii).time, unwrap(res(ii).Z(:, whichvars)), tt), ...
		'movmean', 100);
% 	polarplot(ax(ii), data, tt, '-', 'linewidth', 2);
	plot3(ax(ii), cos(data), sin(data), tt); hold on
	ax(ii).ColorOrderIndex = 1;
	plot3(ax(ii), cos(res(ii).Z(:, whichvars)), ...
		sin(res(ii).Z(:, whichvars)), res(ii).time, '.', 'markersize', 10); hold off
% 	title([pat ' ' res(ii).name]);
	axis tight;
% 	hold on;
% 	set(gca, 'ColorOrderIndex', 1)
	ax(ii + nF) = polaraxes();
	subplot(2, nF, ii + nF, ax(ii + nF));
	polarplot(ax(ii + nF), res(ii).Z(:, whichvars), res(ii).time, '.', ...
		'markersize', 10);
	hold off;
	axis tight
end
legend(fields(whichvars), 'position', [.9 .55 0 0])
rlim = arrayfun(@(a) a.RLim(2), ax);
for ii = 1:numel(ax)
	ax(ii).RLim = [0 max(rlim)];
	ax(ii).ThetaTickLabel = [];
	ax(ii).RTickLabel = [];
end

