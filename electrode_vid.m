function [v] = electrode_vid(data, position, Time, endTime, FrameRate, outfile, visible)

CREATEVID = true;
if ~exist('outfile', 'var') || isempty(outfile)
	outfile = 'mea';
end

if ~exist('FrameRate', 'var') || isempty(FrameRate)
	FrameRate = 1 / (Time(2) - Time(1));
end

if ~exist('visible', 'var') || isempty(visible)
    visible = false;
end

if ~iscell(data)
	data = {data};
end

if ~iscell(position)
	temp = cell(size(data)); 
	[temp{:}] = deal(position); 
	position = temp;
end
numPlots = numel(data);
rows = floor(sqrt(numPlots));
cols = ceil(numPlots / rows);

if numel(Time) ~= size(data{1}, 1)
	error('Data and Time must have the same dimensions.')
end


%% Write video
if CREATEVID
	v = VideoWriter(outfile);
	v.FrameRate = FrameRate;
	open(v);
end

map = make_diverging_colormap('cool', .25);
% map = make_diverging_colormap(1 - parula, 1);

figure(1); clf;
if ~visible, set(1, 'visible', 'off'); end
set(1, 'position', [0   0   335 * cols   254 * rows]);
colormap(map)
% 	p1 = subplot(1,5,1:3);

for p = 1:numPlots
	h(p) = subplot(rows, cols, p);
	
	mNP = quantile(data{p}(:), .025);
	% mE = quantile(dataEWsm(:), .05);

	MNP = quantile(data{p}(:), .975);
	% ME = quantile(dataEWsm(:), .95);

	if MNP == mNP
		MNP = mNP + 1;
	end

	
	scatter(position{p}(:, 1), position{p}(:, 2), 150, data{p}(1, :), 's', 'filled')
	set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
	axis image
	xlim([0 max(position{p}(:, 1))+1])
	ylim([0 max(position{p}(:, 2))+1])
	% 	imagesc(squeeze(dataE(t, :, :)), [mE ME]);
	colorbar
	set(gca,'nextplot','replacechildren');
end

for t = 1:length(Time)
	
	if Time(t) < 0
		desc = 'preictal';
	elseif Time(t) > endTime
		desc = 'postictal';
	else
		desc = 'ictal';
	end
	
	for p = 1:numPlots
		scatter(h(p), position{p}(:, 1), position{p}(:, 2), 150, data{p}(t, :), 's', 'filled')
		title(h(p), sprintf('T = %0.2f (%s)', Time(t), desc))
	end
	drawnow();

% 	p2 = subplot(1, 5, 4:5);
% 	scatter(xNP, yNP, 200, dataNPWsm(t, :), 's', 'filled')
% 	set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
% 	axis image
% 	xlim([0 max(xNP)+1])
% 	ylim([0 max(yNP)+1])
% % 	imagesc(squeeze(dataNP(t, :, :)), [mNP MNP]);
% 	colorbar


	if CREATEVID
		figure(1);
		try
			frame = getframe(gcf);
			writeVideo(v, frame);
		catch ME
			try
				writeVideo(v, frame);
			catch ME
				fprintf('%d', t)
			end
		end
	end

	if ~mod(t, 10)
		disp(['t = ', num2str(t)])
	end
end

if CREATEVID
	close(v)
	close(v)
	disp('Video saved and closed')
end
	