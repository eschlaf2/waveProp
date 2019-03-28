function [] = png2avi(workingDir)
% writes all PNG files in a given directory to an AVI video

imageNames = dir(fullfile(workingDir, '*.png'));
imageNames = {imageNames.name}';

% name = strsplit(workingDir, filesep);
% name = name{end};

outputVideo = VideoWriter([workingDir '.avi']);
% outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir, imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)
