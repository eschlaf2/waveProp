function [data] = data2grid(talldata, position)
% reshapes mea data from form time x electrode to shape Xelectrode x
% Yelectrode x time

nt = size(talldata, 1);  % Number of time points
dims = max(position);  % electrodes in each dimension
[tt, p1] = ndgrid(1:nt, position(:, 1));  % time and x-position grids
[~, p2] = ndgrid(1:nt, position(:, 2));  % ... and y-position grid

data = nan([dims nt]);  % initialize

inds = sub2ind(size(data), p1(:), p2(:), tt(:));  % get addresses
data(inds) = talldata(:);  % assign input to output