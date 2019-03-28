function [mov] = importAvi(filename, framerate)
% Imports an AVI file to a movie structure. Use imshow(mov(x).cdata) to
% show frame x or movie(gcf, mov, n, framerate) to loop through the video n
% times.

vid = VideoReader(filename);
ii = 1; 
while hasFrame(vid)
	mov(ii) = im2frame(readFrame(vid)); 
	ii = ii + 1; 
end


if ~exist('framerate', 'var') || isempty(framerate)
	framerate = vid.FrameRate;
end
	
figure(); movie(gcf, mov, 1, framerate)