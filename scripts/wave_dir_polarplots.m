% pat = 'c7';
if ~exist('files', 'var'); files = dir([pat '_Seizure*10_wave_prop_1.mat']); end
if ~exist('metrics', 'var'); metrics = {'maxdescent', 'events', 'delays_T01_fband1_50'}; end
if ~exist('plotnum', 'var'); plotnum = 1; end
if ~exist('sig', 'var'); sig = 5e-2; end

switch plotnum
    case 1
        nF = numel(files);

        clear res;
        res(nF) = struct(...
            'name', [], ...
            'data', [], ...
            'Z', [], ...
            'time', [], ...
            'p', [], ...
            'Vx', [], ...
            'Vy', []);

        figure(1); clf; fullwidth(true);
        ax = gobjects(2*nF, 1);
        for ii = 1:nF
            res(ii).name = strrep(files(ii).name(strfind(files(ii).name, 'Seizure')+(7:8)), '_', '');
			res(ii).data = load(fullfile(files(ii).folder, files(ii).name), metrics{:});
% 			for m = metrics
% 				res(ii).data.(m{:}) = temp.(m{:});
% 			end
            fields = fieldnames(res(ii).data);
%             whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
			alltimes = cellfun(@(f) ...
				res(ii).data.(f).computeTimes(:), ...
				fields, 'uni', 0);
			alltimes = unique(cat(1, alltimes{:}));
            res(ii).time = alltimes / 1e3;
            ax(ii) = subplot(2, nF, ii, polaraxes());
            ax(ii + nF) = subplot(2, nF, ii + nF, polaraxes());
            [res(ii), ax(ii), ax(ii + nF)] = plot_wave_polar(res(ii), sig, ax(ii), ax(ii + nF));

        end
        legend(strrep(fields, '_', ' '), 'position', [.9 .55 0 0])
        rlim = arrayfun(@(a) a.RLim(2), ax);
        for kk = 1:numel(ax)
            ax(kk).RLim = [0 max(rlim)];
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

%% quantify distributions		
	case 6 
		nF = numel(files);
		metricpairs = nchoosek(metrics, 2);
		nM = size(metricpairs, 1);
		
% 		h = figure(); 
		nrows = nM * nF;
		filename = cell(nrows, 1);
		whichpair = zeros(nrows, 1, 'uint16');
		dZ = cell(nrows, 1);
		[m1, R, theta] = deal(zeros(nrows, 1));
		
		idx = 0;
		for f = 1:nF % for each file
			
			
			[~, name, ~] = fileparts(files(f).name);
			data = load(name);
			
			for m = metricpairs'  % and each pair of metrics
				idx = idx + 1;
				tt = data.(m{1}).computeTimes;
				d1 = data.(m{1}).Z(:);
				d2 = interp1(data.(m{2}).computeTimes, data.(m{2}).Z, tt(:), ...
					'nearest');
				
				dZ{idx} = exp(1j * (d2 - d1));
				m1(idx) = mean(dZ{idx}, 'omitnan');
				R(idx) = abs(m1(idx));  % recall circular variance is 1 - R
				theta(idx) = angle(m1(idx));
				whichpair(idx) = mod(idx-1, nM) + 1;
				filename{idx} = name;
				
			end
			
		end
		

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
            title(ax2, strrep(sprintf('%s - %s', metrics{d1}, metrics{d2}), '_', ' '))
            axis tight
            ax2.ThetaTickLabel = [];
        %     set(gca, 'thetaticklabels', [])
        end
        ttl(sprintf('%s', pat));


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
		K = gausswin(k) * gausswin(k)';  
		% K = ones(k) / k^2;
		im = zeros(dim);
		im(inds) = 1;
		imagesc(ax, conv2(im, K, 'same'))
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

