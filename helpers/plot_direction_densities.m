function plot_direction_densities(res, outer)
% <outer> can be 'files' or 'metrics'

if ~exist('outer', 'var'), outer = 'metrics'; end
metrics = fieldnames(res(1).data);
nF = numel(res);

if strcmpi(outer, 'files')
	OUTER = cellfun(@(f) f(strfind(f, ' ')+1:end), {res.name}, 'uni', 0);
	INNER = rename_metrics(metrics);
else
	OUTER = rename_metrics(metrics);
	INNER = cellfun(@(f) f(strfind(f, ' ')+1:end), {res.name}, 'uni', 0);
end


edges = linspace(-pi, pi, 65);
halfwin = 5;  % (s)
[N, time, Z] = deal(cell(nF, 1));

for iF = 1:nF
	tt = res(iF).time;
	time{iF} = tt(1):.1:tt(end);
	% Rotate Z so that the distributions are centered and it is easy to see
	% deviations
% 			rotateby = angle(sum(exp(1j*res(iF).Z), 'omitnan'));
	rotateby = 0;
	Z{iF} = angle(exp(1j * (res(iF).Z - rotateby)));
	[~, nM] = size(Z{iF});
	N{iF} = nan(numel(time{iF}), numel(edges) - 1, nM);

	for iM = 1:nM
		for tidx = 1:numel(time{iF})
			mask = (tt >= (time{iF}(tidx) - halfwin)) & (tt <= (time{iF}(tidx) + halfwin));
			N{iF}(tidx, :, iM) = histcounts(Z{iF}(mask, iM), edges);
		end

	end
end

nO = numel(OUTER);
nI = numel(INNER);
h = gobjects(nO, 1);

for iO = 1:nO
	h(iO) = figure('Units', 'normalized', 'Position', [0 0 .5 1]);
end


for iO = 1:nO

	patname = strsplit(res(iF).name);

	figure(h(iO));
	annotation('textbox', ...
		'String', [patname{1} ' ' strrep(OUTER{iO}, '_', ' ')], ...
		'Position', [0 0 1 1], ...
		'FitBoxToText', true, ...
		'LineStyle', 'none');

	for ii = 1:nI
		h(iO); subplot(nI, 4, ii * 4 - (3:-1:1)), 
		if strcmpi(outer, 'files')
			iF = iO; 
			iM = ii;
		else
			iF = ii; 
			iM = iO;
		end
		win = gausswin(round(2*halfwin / diff(time{iF}(1:2)))) * ...
			gausswin(round(pi / 2 / diff(edges(1:2))))';
		Zdens = conv2(repmat(N{iF}(:, :, iM), 1, 3), win, 'same');
		Zdens = Zdens(:, size(N{iF}, 2):2*size(N{iF}, 2));
		Zdens = Zdens ./ max(sum(Zdens, 2), mean(sum(Zdens, 2))) * 100;
		emilys_pcolor(time{iF}, edges, Zdens, 'clims', [0 5]); %, ...
% 			'clims', [0 20]); 
		hold on;
		try
			plot(res(iF).time, (angle((exp(1j * Z{iF}(:, iM))))), 'c.');
		catch ME
			if ~strcmp(ME.identifier, 'MATLAB:griddedInterpolant:DegenerateGridErrId')
				rethrow(ME)
			end
			disp(ME)
		end
		hold off;
		ytks = -pi:pi/2:pi;
		ytklbl = cell(size(ytks));
		ytklbl{2} = '-\pi/2'; ytklbl{4} = '\pi/2';
		for y=ytks, yline(y); end  % grid doesn't show over pcolor

		yticks(ytks);
		yticklabels(ytklbl);
		xlim([min(cat(2, time{:})), max(cat(2, time{:}))]);
		title(strrep(INNER{ii}, '_', ' '))

		subplot(nI, 4, ii * 4);
		plot(max(Zdens), edges, 'linewidth', 2); axis tight;
		yticks(ytks);
		yticklabels([]);
		xlim([0 inf]);
		grid on;
	end


	drawnow;
	print(sprintf('figs/dens_%s_%s', patname{1}, checkname(OUTER{iO})), '-dpng');
end
