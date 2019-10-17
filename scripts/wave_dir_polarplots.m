% pat = 'c7';
% Requires Circular Statistics Toolbox (plotnum case 6)

if ~exist('seizure', 'var'); seizure = '*'; else, seizure = num2str(seizure); end
if ~exist('files', 'var'); files = dir([pat '_Seizure' seizure '*10_wave_prop_1.mat']); end
if ~exist('metrics', 'var')
	metrics = {...
		'delays_T0p2_fband0_50', ...
		'maxdescent', ...
		'events', ...
		'delays_T01_fband1_50', ...
		'delays_T10_fband1_13'}; 
end
if ~exist('plotnum', 'var'); plotnum = 0; end
if ~exist('sig', 'var'); sig = 5e-2; end
nF = numel(files);

%% Make res


switch plotnum
	case 0
	% Make res
		clear res;
        res(nF) = struct(...
            'name', [], ...
            'data', [], ...
            'Z', [], ...
            'time', [], ...
            'p', [], ...
            'Vx', [], ...
            'Vy', []);

        for ii = 1:nF
			finfo = strsplit(files(ii).name, {'_', '.'});
            res(ii).name = sprintf('%s %s', finfo{1}, finfo{2}(8:end));
			res(ii).data = load(fullfile(files(ii).folder, files(ii).name), metrics{:});
            fields = fieldnames(res(ii).data);
			alltimes = cellfun(@(f) ...
				res(ii).data.(f).computeTimes(:), ...
				fields, 'uni', 0);
			alltimes = cat(1, alltimes{:});
            res(ii).time = unique(alltimes / 1e3, 'stable');
			
			
			time = res(ii).time;  % - min(res.time);
			[res(ii).Z, res(ii).Vx, res(ii).Vy] = ...
				deal(zeros(length(time), length(fields)));
			for jj = 1:numel(fields)
				for f = 'Zpxy'

					switch f
						case {'Z', 'p'}
							data = interp1(...
								res(ii).data.(fields{jj}).computeTimes / 1e3, ...
								res(ii).data.(fields{jj}).(f), time, 'nearest');
						case 'x'
							data = interp1(...
								res(ii).data.(fields{jj}).computeTimes / 1e3, ...
								res(ii).data.(fields{jj}).V(1, :), time, 'nearest');
							f = 'Vx';
						case 'y'
							data = interp1(...
								res(ii).data.(fields{jj}).computeTimes / 1e3, ...
								res(ii).data.(fields{jj}).V(2, :), time, 'nearest');
							f = 'Vy';
					end
					res(ii).(f)(:, jj) = data;
					% Remove values where fit is not significant
					mask = res(ii).data.(fields{jj}).p < sig;
					res(ii).(f)(~mask, jj) = nan;
					% Remove values where slope is zero in both directions
					res(ii).(f)(all(abs(res(ii).data.(fields{jj}).beta(1:2, :)) < eps), jj) = nan;
				end
			end
        end
        
		
    case 1

        figure(1); clf; fullwidth(true);
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
        ttl(pat)
		
%% Density plots
	case 7
		% Plots histogram (normalized to percentage) of directions within 5 seconds
% of each time point at .1 second intervals. Each seizure and metric

edges = linspace(-pi, pi, 65);
halfwin = 5;  % (s)
N = cell(nF, 1);
for ii = 1:nF
	tt = res(ii).time;
	time = tt(1):.1:tt(end);
	% Rotate Z so that the distributions are centered and it is easy to see
	% deviations
	rotateby = angle(sum(exp(1j*res(ii).Z), 'omitnan'));
	Z = angle(exp(1j * res(ii).Z) .* exp(-1j * rotateby));
	[~, nM] = size(Z);
	N{ii} = nan(numel(time), numel(edges) - 1, nM);
	figure('Units', 'normalized', 'Position', [0 0 .5 1])
	win = gausswin(round(2*halfwin / diff(time(1:2)))) * ...
		gausswin(round(pi / 2 / diff(edges(1:2))))';
	for m = 1:nM
		for tidx = 1:numel(time)
			mask = (tt >= (time(tidx) - halfwin)) & (tt <= (time(tidx) + halfwin));
			N{ii}(tidx, :, m) = histcounts(Z(mask, m), edges);
		end
		subplot(nM, 1, m), 
		Zdens = conv2(N{ii}(:, :, m), win, 'same');
		Zdens = Zdens ./ max(sum(Zdens, 2), mean(sum(Zdens, 2))) * 100;
		emilys_pcolor(time, edges(2:end), Zdens, 'clims', [0 5]); %, ...
% 			'clims', [0 20]); 
		hold on;
		try
		plot(tt, (angle(smoothdata(exp(1j * Z(:, m)), 'omitnan', ...
			'SamplePoints', tt))), 'c.');
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
		title(rename_metrics(metrics{m}))
	end
	annotation('textbox', ...
		'String', res(ii).name, ...
		'Position', [0 0 1 1], ...
		'FitBoxToText', true, ...
		'LineStyle', 'none');
	drawnow;
	print(sprintf('figs/dens_%s', checkname(res(ii).name)), '-dpng');
end


%% quantify distributions		
	case 6 

		metricpairs = nchoosek(metrics, 2);
		nM = size(metricpairs, 1);
		
		nrows = nM * nF;
		filename = cell(nrows, 1);
		whichpair = zeros(nrows, 1, 'uint16');
		dZ = cell(nrows, 1);
		[m1, R, theta, kappa, conf, sigma, N] = ...
            deal(nan(nrows, 1));
		
		idx = 0;
		for f = 1:nF % for each file
			
			
			[~, name, ~] = fileparts(files(f).name);
			data = load(name);
			
			for m = metricpairs'  % and each pair of metrics
				idx = idx + 1;
				
				whichpair(idx) = mod(idx-1, nM) + 1;
				filename{idx} = name;
				
				tt = data.(m{1}).computeTimes;
				d1 = data.(m{1}).Z(:);
                d1(data.(m{1}).p(:) >= sig) = nan;
                d2 = data.(m{2}).Z(:);
                d2(data.(m{2}).p(:) >= sig) = nan;
				d2 = interp1(data.(m{2}).computeTimes, d2, tt(:), ...
					'nearest');
				dd = d2 - d1;
				dZ{idx} = exp(1j * (dd));
				finite = isfinite(dd);
                N(idx) = sum(finite);
				if N(idx) == 0, continue, end
                
				m1(idx) = mean(dZ{idx}, 'omitnan');
                R(idx) = abs(m1(idx));
                kappa(idx) = circ_kappa(dd(finite));
                conf(idx) = circ_confmean(dd(finite));
                sigma(idx) = circ_std(dd(finite));
				theta(idx) = angle(m1(idx));
				
			end
			
		end
        stats = table(filename, whichpair, R, theta, kappa, conf, sigma, N, ...
            m1, dZ);
		

%% Comparison of metrics

    case 2
        alltime = cell(numel(metrics), 1);
        % Make a matrix of all metrics from all seizures
        for seiz = 1:nF  % for each seizure
            timetemp = res(seiz).time;
            alltime{seiz} = (timetemp - min(timetemp)) / range(timetemp);
            fields = fieldnames(res(seiz).data);
            temp = res(seiz).Z(:, strcmpi(metrics{1}, fields));  % Get the first metric
            for m = 2:numel(metrics)  % concatenate the remaining metrics as columns
                temp = cat(2, temp, res(seiz).Z(:, strcmpi(metrics{m}, fields)));
            end
            if seiz == 1  % If this is the first seizure
                alldata = temp;  % ... initialize alldata
            else  % otherwise
                alldata = cat(1, alldata, temp);  % ... concat data to rows
            end
        end
        alltime = cat(1, alltime{:});  % make a long column of all time points 
        % ... in time relative to seizure progression


        metricpairs = nchoosek(1:numel(metrics), 2);  % Get each pair of metrics
        nP = size(metricpairs, 1);
        r = ceil(sqrt(nP)); c = ceil(nP / r);  % Create an appropriate number of subplots
        figure(2); clf; fullwidth(r > 1);
        for dd = 1:nP  % For each pair
            ax2 = polaraxes();  % Create a polar axis,
            subplot(r, c, dd, ax2)  % ... place it,
            d1 = metricpairs(dd, 1); d2 = metricpairs(dd, 2);  % ... assign convenience vars,
            data = alldata(:, d1) - alldata(:, d2);  % ... compute the difference between angles from each metric,
            edges = linspace(0, 2*pi, 41);  % ... create 40 bins between 0 and 2pi,
            bc = histcounts(mod(data, 2*pi), edges);  % ... get the counts per bin,
            bc = bc/max(bc);  % ... scale to match time (relative), 
            if ~any(isfinite(bc)); continue; end
            polarhistogram(ax2, 'BinEdges', edges, 'BinCounts', bc); hold on;  % ... draw the histogram,
            polarplot(ax2, data, alltime, '.'); hold off;  % ... add the differences at each time

            % Prettify
            title(ax2, strrep(sprintf('%s - %s', ...
				rename_metrics(metrics{d1}), ...
				rename_metrics(metrics{d2})), '_', ' '))
            axis tight
            ax2.ThetaTickLabel = [];
        %     set(gca, 'thetaticklabels', [])
		end
		
		annotation('textbox', ...
			'String', pat, ...
			'Position', [0 0 1 1], ...
			'FitBoxToText', true, ...
			'LineStyle', 'none');
		drawnow;
%         ttl(sprintf('%s', pat));


%% Blurred im
HIDE = true;

if ~HIDE
	dim = 100;
	figure(3); fullwidth()
	scale = @(x, M) round((x - min(x)) / range(x) * (M - 1) + 1);
	for file = 1:nF
		fields = fieldnames(res(file).data);
		whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
		ax = subplot(1, nF, file);
		time = res(file).time(:); time = time - min(time);
		xdata = time .* cos(res(file).Z(:, whichfields)); 
		ydata = time .* sin(res(file).Z(:, whichfields)); 
		valid = ~isnan(xdata);
		cdata = ones(size(ydata));
		inds = sub2ind([dim, dim], ...
			scale(xdata(valid), dim), scale(ydata(valid), dim));
		k = 10;
		kappa = gausswin(k) * gausswin(k)';  
		% K = ones(k) / k^2;
		im = zeros(dim);
		im(inds) = 1;
		imagesc(ax, conv2(im, kappa, 'same'))
		axis square
		axis xy
	end
end
%% Plot p-values
    case 4
    figure(4); fullwidth();

    seizure = 1;
    thresh = 5e-10;

    whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
    p = res(seizure).p(:, whichfields);
    p(p > thresh) = nan;
    stem(res(seizure).time, -log(p), 'filled', 'linewidth', 2)

    legend(metrics)
    title('p-values')

%%
    case 5
        figure(5); fullwidth(true)
        seizure = 3;
        whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));

        a1 = subplot(2,1,1);
        stem(res(seizure).time, res(seizure).Z(:, whichfields), 'filled');
        legend(metrics)

        a2 = subplot(2,1,2);
        d1 = -diff(res(seizure).Z(:, whichfields)');
        d1(d1 > pi) = d1(d1 > pi) - 2*pi;
        d1(d1 < -pi) = d1(d1 < -pi) + 2*pi;
        stem(res(seizure).time, d1, 'filled');
        hold on
        plot(res(seizure).time, smoothdata(d1, 'movmean', 10))
        hold off;
        legend('difference')

        linkaxes([a1, a2], 'xy')

%        print(2, [pat '_polarhist'], '-dpng');
%        print(5, [pat '_Seizure' res(seizure).name '_diff'], '-dpng')


end

