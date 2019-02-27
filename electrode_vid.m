function [v] = electrode_vid(temp, X, Y, Time, endTime, FrameRate, outfile, visible)

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

if numel(Time) ~= size(temp, 1)
	error('Data and Time must have the same dimensions.')
end
%% Reshape data
% tmax = length(Time);
% dataR = nan(tmax, max(X), max(Y));
% 
% [tN, xN, yN] = size(dataR);
% 
% 
% for t = 1:tmax
% 	inds = sub2ind([tN, xN, yN], t * ones(size(X)), X, Y);
% 	dataR(inds) = temp(t, :);
% % 	dataR(t + (X - 1) * tN + (Y - 1) * tN * xN) = ...
% % 		temp(t, :);
% end

%% Write video
if CREATEVID
	v = VideoWriter(outfile);
	v.FrameRate = FrameRate;
	open(v);
end

mNP = quantile(temp(:), .05);
% mE = quantile(dataEWsm(:), .05);

MNP = quantile(temp(:), .95);
% ME = quantile(dataEWsm(:), .95);

if MNP == mNP
	MNP = mNP + 1;
end

map = make_diverging_colormap('cool', 1);

figure(1); clf;
if ~visible, set(1, 'visible', 'off'); end
set(1, 'position', [556   422   335   254]);
colormap(map)
% 	p1 = subplot(1,5,1:3);
scatter(X, Y, 150, temp(1, :), 's', 'filled')
set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
axis image
xlim([0 max(X)+1])
ylim([0 max(Y)+1])
% 	imagesc(squeeze(dataE(t, :, :)), [mE ME]);
colorbar
set(gca,'nextplot','replacechildren');

for t = 1:length(Time)
	
	if Time(t) < 0
		desc = 'preictal';
	elseif Time(t) > endTime
		desc = 'postictal';
	else
		desc = 'ictal';
	end
	
	
	scatter(X, Y, 150, temp(t, :), 's', 'filled')
	title(sprintf('T = %0.2f (%s)', Time(t), desc))
% 	drawnow();

% 	p2 = subplot(1, 5, 4:5);
% 	scatter(xNP, yNP, 200, dataNPWsm(t, :), 's', 'filled')
% 	set(gca, 'clim', [mNP MNP], 'Color', .15*[1 1 1]);
% 	axis image
% 	xlim([0 max(xNP)+1])
% 	ylim([0 max(yNP)+1])
% % 	imagesc(squeeze(dataNP(t, :, :)), [mNP MNP]);
% 	colorbar


	if CREATEVID
		try
			frame = getframe(gcf);
			writeVideo(v,frame);
		catch ME
			try
				writeVideo(v, frame);
			catch ME
				fprintf('%d', t)
			end
% 			close(v);
% 			disp(t);
% 			rethrow(ME);
		end
	end
% 		drawnow()
% 	end

	if ~mod(t, 10)
		disp(['t = ', num2str(t)])
	end
end

if CREATEVID
	close(v)
	close(v)
	disp('Video saved and closed')
end
	