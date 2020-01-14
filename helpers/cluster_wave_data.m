function T = cluster_wave_data(data, position)
% Improve wave fitting - e.g. if there are two waves moving in different
% directions, fitting a plane wave can give a bad result. 
% Inputs:
%	data: nx1 time of occurrence
%   position: nx2 electrode locations
 
data = [data(:) position];

%%% Underlying algorithm
% Y = pdist(data, 'seuclidean');  
% Z = linkage(Y, 'single');
% T = cluster(Z, 'cutoff', .5, 'criterion', 'distance');

T = clusterdata(data, ...
	'distance', 'seuclidean', ... % standardized euclidean distance
	'cutoff', .7, 'criterion', 'distance');  % Separate large gaps (this distance is somewhat arbitrary right now...)

%%% For copy-paste on SCC
% T = clusterdata(data, 'distance', 'seuclidean','cutoff', .5, 'criterion', 'distance'); 
% figure(5); scatter3(position(:, 1), position(:, 2), wave_fit.data{w}, 30, T, 'filled')
