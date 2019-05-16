function ax = plot_window(window, position, ax, v)

if ~exist('ax', 'var') || isempty(ax)
	figure(); ax = axes();
end
title_string = ax.Title.String;

writevid = true;
if ~exist('v', 'var') || isempty(v)
	writevid = false;
end

[nT, ~] = size(window);  % First dimension should be time in window
position = (position - min(position) + 1) ./ min(diff(unique(position)));  % Reset position to be integers 1 .. n
addy = sub2ind(max(position), position(:, 1), position(:, 2));  % convert to indices
clims = [quantile(window(:), .025) quantile(window(:), .975)];  % Set color limits based on entire window
zlims = [min(min(window(:)), 0), max(window(:))];  % Set z limits based on entire wave
for jj = 1:nT  % show data at each time step
	tempV = nan(max(position));  % make a clean array
	tempV(addy) = window(jj, :)';  % fill in the data
	surf(ax, tempV, 'FaceAlpha', .8);
	zlim(zlims);
	hold(ax, 'on');
	imagesc(ax, tempV, clims); colorbar(ax);  % display
	hold(ax, 'off');
	xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'lfp');
	title(title_string);
	if writevid
% 		drawnow()
		f = getframe(ax.Parent);
		writeVideo(v, f)
	else
		drawnow()
	end
end
