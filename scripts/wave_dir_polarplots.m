% pat = 'c7';
% Requires Circular Statistics Toolbox (plotnum case 6)
% fid = fopen('seizures2.txt');
% A = textscan(fid, '%s %d'); 
% pats = A{1}; seizures = A{2};
% for ii = numel(pats):-1:1
% 	files(ii) = dir(sprintf('%s_Seizure%d_Neuroport_10_10_wave_prop.mat', pats{ii}, seizures(ii)));
% end
% fclose(fid);

% if ~exist('seizure', 'var'); seizure = '*'; elseif isnumeric(seizure), seizure = num2str(seizure); end
% % if ~exist('files', 'var'); files = dir([pat '_Seizure' seizure '_Neuroport_10_10_wave_prop.mat']); end
% if ~exist('metrics', 'var')
% 	metrics = {...
% 		'maxdescent', ...
% 		'events', ...
% 		'delays_T10_fband1_13', ...
% 		'delays_T01_fband1_13'}; %, ...
% % 		'delays_T10_fband1_50', ...
% % 		'delays_T01_fband1_50'}; 
% end
% if ~exist('plotnum', 'var'); plotnum = 0; end
% if ~exist('sig', 'var'); sig = 5e-2; end
% nF = numel(files);

%% Make res


switch plotnum
	case 0
	% Make res
		res = compile_wave_prop();
		
    case 1

        figure(1); clf; fullwidth(true);
		plot_circular_trajectories(res);
		

%% Density plots
	case 7
		% Plots histogram (normalized to percentage) of directions within 5 seconds
		% of each time point at .1 second intervals. Each seizure and metric
		
		outer = 'metrics';  % 'metrics' or 'files'
		
		plot_direction_densities(res, outer)
		

%% quantify distributions		
	case 6 

		if ~exist('res', 'var'); res = compile_wave_prop; end
		if ~exist('metrics', 'var'), metrics = fieldnames(res(1).data); end
		if ~exist('sig', 'var'), sig = 5e-2; end
		nF = length(res);
		metricpairs = nchoosek(metrics, 2);
		nM = size(metricpairs, 1);
		
		nrows = nM * nF;
		name = cell(nrows, 1);
		whichpair = zeros(nrows, 1, 'uint16');
		dZ = cell(nrows, 1);
		[m1, R, theta, kappa, conf, sigma, N] = ...
            deal(nan(nrows, 1));
		
		idx = 0;
		for f = 1:nF % for each file
			
			data = res(f).data;
			for m = metricpairs'  % and each pair of metrics
				idx = idx + 1;
				
				whichpair(idx) = mod(idx-1, nM) + 1;
				name{idx} = strsplit(res(f).name);
				
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
        stats = table(name, whichpair, R, theta, kappa, conf, sigma, N, ...
            m1, dZ);
		save('direction_stats', 'stats', 'metricpairs');
		

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

    s = 1;

    whichfields = find(sum(cell2mat(cellfun(@(f) strcmpi(f, fields), metrics, 'uni', 0)), 2));
    p = res(s).p(:, whichfields);
    p(p > sig) = nan;
%     plot(res(s).time, -log(p), '.', 'markersize', 20)
	bar(res(s).time, -log(p), 'stacked')

    legend(rename_metrics(metrics))
    title([pat ' p-values'])
	ylabel('-log(p)')
	xlabel('Time (s)')

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

