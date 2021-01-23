% toi = [50 60];
% fname = 'MG49_Seizure43_Neuroport_10_10.mat';
% skipfactor = 30;
% clims = [];
% bands = [[2; 2.2] [12; 12.2] [22; 22.2]];
% T = 5; step = .1;
% framerate = 100;
% whiten = 0;

% mea = MEA(fname);
% mea.whiten;
disp(mea.Name);
pat = mea.patient; seizure = str2double(mea.seizure);
if ~isempty(clims), climfun = @(data, ii) quantile(single(data{ii}(:)), clims); end

% mea = exclude_channels(mea);

%%
position = mea.Position; 
mea.SamplingRate = framerate;
time = mea.Time(); timeinds = find(time > toi(1) & time < toi(2));
doi = smoothdata(single(mea.Data(timeinds, :)));

%%
mea.params.fr_window = 20;
fr = mea.firing_rate(timeinds, :);


%%
clear data
N = length(doi);
% doiR = nan([max(position) N]);
% [tt, p1] = ndgrid(1:N, position(:, 1));
% [~, p2] = ndgrid(1:N, position(:, 2));
% addy = sub2ind(size(doiR), p1(:), p2(:), tt(:));
% doiR(addy) = doi(:);
data{1} = doi;
ttl{1} = 'raw';

data{2} = fr;
ttl{2} = 'Firing rate';

    
ii = 3;
for passband = bands
    [dat_temp, pb] = mea.filter(mea.Raw, mea.SRO, passband); 
    dat_temp = single(resample(double(dat_temp), 1, mea.skipfactor));
    data{ii} = dat_temp(timeinds, :);
    ttl{ii} = sprintf('Filtered to [%.1f %.1f]', pb);
    ii = ii + 1;
end

% t0 = time(1); tf = t0 + T;
% for passband = bands
% 	band = passband;
% 	if band(1) == 0, band(1) = []; end
% 	band = min(band, samplerate / 2 - 1);
%     b = fir1(150, 2 * band / samplerate);  % convert passband to pi-rad/sample
%     tempdat = zeros(size(doi));
%     while tf <= time(end)
%         inds = (time >= t0) & (time <= tf);
%         tempf = single(filtfilt(b, 1, double(doi(inds, :))));
%         tempdat(inds, :) = (tempdat(inds, :) + tempf) / 2;
%         t0 = t0 + step; tf = t0 + T;
%     end
%     
%     data{ii} = single(filtfilt(b, 1, double(doi)));
%     ttl{ii} = sprintf('Filtered to [%.1f %.1f]', passband);
%     ii = ii + 1;
% end

%%
close all
numplots = length(data);
r = floor(sqrt(numplots)); c = ceil(numplots / r);
h = figure(1); clf; set(1, 'Position', [0 0 300 * c 225 * r]); colormap(1-gray)

ax = gobjects(numplots, 1);
for ii = 1:numplots
	ax(ii) = subplot(r, c, ii); 
    if ~isempty(clims), ax(ii).CLim = climfun(data, ii); end
	title(ax(ii), ttl{ii})
	ax(ii).XLim = [0 max(position(:, 1)) + 1];
	ax(ii).YLim = [0 max(position(:, 2)) + 1];
    ax(ii).XTickLabel = [];
    ax(ii).YTickLabel = [];
	ax(ii).NextPlot = 'replacechildren';
	if ~isempty(clims), colorbar(ax(ii)); end
    scatter(ax(ii), ...
        position(:, 1), position(:, 2), 225, data{ii}(1, :), ...
        'filled', 'square')

end
% sgtitle(h, sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, toi(1)));
th = annotation(h, 'textbox', [0 1 0 0], 'String', sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, toi(1)), 'FitBoxToText', 'on', 'LineStyle', 'none');
clear mov
mov(N) = getframe(h);

%%
inds = 1:N;
for ii = inds  
	if ~mod(ii, 100), fprintf('Writing frame %d/%d\n', ii, inds(end)), end
	for p = 1:numplots
        ax(p).Children.CData = data{p}(ii, :);
	end
	th.String = sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, time(ii));
	drawnow
	mov(ii) = getframe(h);
end

if all(isinf(toi))
    t_str = '';
else
    t_str = strrep(sprintf('_time_%03.0f_%03.0f', toi), '-', 'M');
end
if isempty(clims)
    filename = sprintf('%s%s_climsAuto', mea.Name, t_str);
else
    filename = sprintf('%s%s_clims%03.0f_%03.0f', ...
        mea.Name, t_str, clims(1)*100, clims(2)*100);
end

if whiten
    v = VideoWriter([filename '_whitened']);
else
    v = VideoWriter(filename);
end
disp(['Saving ' v.Filename ' ...'])
v.FrameRate = 30;
open(v);
writeVideo(v, mov(inds));
close(v);
disp('Done.')
