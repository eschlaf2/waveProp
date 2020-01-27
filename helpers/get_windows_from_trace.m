function windows = get_windows_from_trace(data, samples, sample_points)
% Extract windows surrounding points in time or sampled at a given
% interval. 
% Inputs:
%	data: t-by-k array
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


if numel(samples) == 1, samples = 0:samples-1; end  % if samples is a length, convert to indices
if numel(sample_points) == 1  % if sample_points is a skipfactor, convert to sample_points
	sample_points = 0:sample_points:length(data) - samples(end); 
end

% [tt, ww] = ndgrid(samples, sample_points);
% inds = tt + ww;
windows = nan(numel(samples), size(data, 2), numel(sample_points));
for ii = 1:numel(sample_points)
	windows(:, :, ii) = data(sample_points(ii) + samples, :);
end
windows = squeeze(windows);


