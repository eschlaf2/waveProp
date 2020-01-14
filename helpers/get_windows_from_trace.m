function windows = get_windows_from_trace(data, samples, sample_points)
% Extract a period surrounding points or time or sampled at a given
% interval. 
% Inputs:
%	data: nx1 vector or the data
%	samples: (vector) indices surrounding sample sample_points to return OR
%			 (scalar) number of points to collect starting at sample_points
%					  E.g. the following are the same:
%							samples = N;
%							samples = 0:N-1;
%	sample_points: (vector) indices at which to sample the data
%				   (scalar) skipfactor - sample every skipfactor points
% 					  E.g. the following are the same:
%							sample_points = N;
%							sample_points = 1:N:length(data)-samples(end)


data = data(:);
if numel(samples) == 1, samples = 0:samples-1; end  % if samples is a length, convert to indices
if numel(sample_points) == 1  % if sample_points is a skipfactor, convert to sample_points
	sample_points = 0:sample_points:length(data) - samples(end); 
end

[tt, ww] = ndgrid(samples, 1:sample_points:length(data) - samples(end));
inds = tt + ww;
windows = reshape(data(inds), size(tt));
