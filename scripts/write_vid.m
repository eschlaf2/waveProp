% toi = [50 60];
% fname = 'MG49_Seizure43_Neuroport_10_10.mat';
% skipfactor = 30;
% clims = [];
% bands = [[2; 2.2] [12; 12.2] [22; 22.2]];
% T = 5; step = .1;

mea = load(fname);
disp(mea);
[~, name, ~] = fileparts(fname);
finfo = strsplit(name, {'_', filesep});
pat = finfo{1}; seizure = str2double(finfo{2}(8:end));
if ~isempty(clims), climfun = @(data, ii) quantile(single(data{ii}(:)), clims); end

% mea = exclude_channels(mea);

%%
position = mea.Position; position(mea.BadChannels, :) = [];
time = mea.Time(); timeinds = find(time > toi(1) & time < toi(2));
time = time(timeinds(1:skipfactor:end));
doi = smoothdata(single(mea.Data(timeinds(1:skipfactor:end), :)));
samplerate = mea.SamplingRate / skipfactor;
doi(:, mea.BadChannels) = [];

%%
clear data
N = length(doi);
doiR = nan([max(position) N]);
[tt, p1] = ndgrid(1:N, position(:, 1));
[~, p2] = ndgrid(1:N, position(:, 2));
addy = sub2ind(size(doiR), p1(:), p2(:), tt(:));
doiR(addy) = doi(:);
data{1} = doi;
ttl{1} = 'raw';

    
ii = 2;
t0 = time(1); tf = t0 + T;
for passband = bands
	band = passband;
	if band(1) == 0, band(1) = []; end
	band = min(band, samplerate / 2 - 1);
    b = fir1(150, 2 * band / samplerate);  % convert passband to pi-rad/sample
    tempdat = zeros(size(doi));
    while tf <= time(end)
        inds = (time >= t0) & (time <= tf);
        tempf = single(filtfilt(b, 1, double(doi(inds, :))));
        tempdat(inds, :) = (tempdat(inds, :) + tempf) / 2;
        t0 = t0 + step; tf = t0 + T;
    end
    data{ii} = single(filtfilt(b, 1, double(doi)));
    ttl{ii} = sprintf('Filtered to [%.1f %.1f]', passband);
    ii = ii + 1;
end

%%
close all
numplots = length(data);
r = floor(sqrt(numplots)); c = ceil(numplots / r);
h = figure(1); clf; set(1, 'Position', [0 0 300 * c 225 * r]); colormap(bone)

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
% inds = find((time >= 42.5) & (time < 42.625));
inds = 1:N;
for ii = inds  % 1:N
	if ~mod(ii, 100), fprintf('Writing frame %d/%d\n', ii, inds(end)), end
	for p = 1:numplots
        ax(p).Children.CData = data{p}(ii, :);
%         ax(p).CLim = [-1400 -250];
%         ax(p).CLim = quantile(double(data{p}(ii, :)), [.01 .15]);
	end
	th.String = sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, time(ii));
	drawnow
	mov(ii) = getframe(h);
end

if isempty(clims)
    filename = strrep(sprintf('%s_time%03.0f_%03.0f_climsAuto', ...
	name, toi(1), toi(2)), '-', 'M');
else
    filename = strrep(sprintf('%s_time%03.0f_%03.0f_clims%03.0f_%03.0f', ...
        name, toi(1), toi(2), clims(1)*100, clims(2)*100), '-', 'M');
end

v = VideoWriter(filename);
disp(['Saving ' v.Filename ' ...'])
v.FrameRate = 30;
open(v);
writeVideo(v, mov(inds));
close(v);
disp('Done.')
