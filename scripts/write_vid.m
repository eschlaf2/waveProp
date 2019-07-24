% toi = [0 10];
% fname = 'MG49/MG49_Seizure43_Neuroport_10_10.mat';
% skipfactor = 1;
% clims = @(data) quantile(single(data(:)), [0 1]);
% bands = [[1; 13] [20; 40]];

mea = load(fname);
disp(mea);
[~, name, ~] = fileparts(fname);
finfo = strsplit(name, {'_', filesep});
pat = finfo{1}; seizure = str2double(finfo{2}(8:end));
% mea = exclude_channels(mea);

%%
position = mea.Position; position(mea.BadChannels, :) = [];
time = mea.Time(); timeinds = time > toi(1) & time < toi(2);
time = downsample(time(timeinds), skipfactor);
doi = downsample(mea.Data(timeinds, :), skipfactor);
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
for passband = bands
    b = fir1(150, 2 * passband / samplerate);  % convert passband to pi-rad/sample
    data{ii} = single(filtfilt(b, 1, double(doi)));
    ttl{ii} = sprintf('Filtered to [%d %d]', passband);
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
	ax(ii).CLim = clims(data{ii});
	title(ax(ii), ttl{ii})
	ax(ii).XLim = [0 max(position(:, 1)) + 1];
	ax(ii).YLim = [0 max(position(:, 2)) + 1];
	ax(ii).NextPlot = 'replacechildren';
	colorbar(ax(ii));
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

v = VideoWriter(sprintf('%s_Seizure%d_Neuroport_10_10_time%03d_%03d', pat, seizure, toi(1), toi(2)));
v.FrameRate = 30;
open(v);
writeVideo(v, mov(inds));
close(v);
