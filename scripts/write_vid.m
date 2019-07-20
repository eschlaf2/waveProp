close all

toi = [45 85];
% toi = [15 25];
pat = 'MG49'; seizure = 43;
% pat = 'c7'; seizure = 1;
skipfactor = 1;

addpath(pat);
mea = load(sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure));
disp(mea);
% mea = exclude_channels(mea);

position = mea.Position; position(mea.BadChannels, :) = [];
time = mea.Time(); timeinds = time > toi(1) & time < toi(2);
time = downsample(time(timeinds), skipfactor);
doi = downsample(mea.Data(timeinds, :), skipfactor);
samplerate = mea.SamplingRate / skipfactor;
doi(:, mea.BadChannels) = [];

%%
clims = quantile(single(doi(:)), single([.01 .15]));
N = length(doi);
doiR = nan([max(position) N]);
[tt, p1] = ndgrid(1:N, position(:, 1));
[~, p2] = ndgrid(1:N, position(:, 2));
addy = sub2ind(size(doiR), p1(:), p2(:), tt(:));
doiR(addy) = doi(:);
data{1} = doi;
ttl{1} = 'raw';

passband = [10 40];
b = fir1(150, 2 * passband / samplerate);  % convert passband to pi-rad/sample
data{2} = single(filtfilt(b, 1, double(doi)));
ttl{2} = sprintf('Filtered to [%d %d]', passband);

% passband = [20 50];
% b = fir1(150, 2 * passband / samplerate);  % convert passband to pi-rad/sample
% data{3} = single(filtfilt(b, 1, double(doi)));
% ttl{3} = sprintf('Filtered to [%d %d]', passband);
% 
% passband = [9 14];
% b = fir1(150, 2 * passband / samplerate);  % convert passband to pi-rad/sample
% data{4} = single(filtfilt(b, 1, double(doi)));
% ttl{4} = sprintf('Filtered to [%d %d]', passband);
% 
% passband = [15 20];
% b = fir1(150, 2 * passband / samplerate);  % convert passband to pi-rad/sample
% data{5} = single(filtfilt(b, 1, double(doi)));
% ttl{5} = sprintf('Filtered to [%d %d]', passband);


numplots = length(data);
r = floor(sqrt(numplots)); c = ceil(numplots / r);
h = figure(1); clf; set(1, 'Position', [0 0 300 * c 225 * r]); colormap(bone)

ax = gobjects(numplots, 1);
for ii = 1:numplots
	ax(ii) = subplot(r, c, ii); 
	ax(ii).CLim = quantile(data{ii}(:), [.01 .99]);
	title(ax(ii), ttl{ii})
	ax(ii).XLim = [0 max(position(:, 1)) + 1];
	ax(ii).YLim = [0 max(position(:, 2)) + 1];
	ax(ii).NextPlot = 'replacechildren';
	colorbar(ax(ii));
% 	h.NextPlot = 'replacechildren';

end
sgtitle(h, sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, toi(1)));
mov(N) = getframe(h);	
for ii = 1:N
	for p = 1:numplots
	scatter(ax(p), ...
		[position(:, 1)], [position(:, 2)], 225, [data{p}(ii, :)], ...
		'filled', 'square')
	end
	h.Children(1).String = sprintf('%s Seizure %d\nT=%.2f\n', pat, seizure, time(ii));
	drawnow
	
	
	mov(ii) = getframe(h);
end

v = VideoWriter(sprintf('%s_Seizure%d_Neuroport_10_10_time%03d_%03d', pat, seizure, toi(1), toi(2)));
v.FrameRate = 30;
open(v);
writeVideo(v, mov);
close(v);
