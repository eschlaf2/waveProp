function [mea, P] = create_mea(Data, varargin)

P = inputParser;
r =@(varargin) addRequired(P, varargin{:});
p =@(varargin) addParameter(P, varargin{:});

% Required
r('Data');

% Optional
p('Time', []);
p('Padding', [10 10]);
p('Duration', []);
p('SamplingRate', Inf);
p('Path', [pwd filesep 'mea.mat']);
p('Position', []);
p('BadChannels', []);
p('Map', []);
p('Name', 'MEA data');
p('firing_rate', []);

parse(P, Data, varargin{:})

mea = P.Results;
tol = 1e-6;

mea = check_data(mea, tol);
mea = check_position(mea, tol);
mea = check_time(mea, tol);

end

function mea = check_time(mea, tol)
% Compare Time and SamplingRate for consistency. 
% Add duration

switch isinf(mea.SamplingRate) + isempty(mea.Time)
	case 2
		error('''SamplingRate'' or ''Time'' must be provided.');
	case 1
		if isinf(mea.SamplingRate)
			t = mea.Time();
			mea.SamplingRate = 1 / diff(t(1:2));
		else
			dt = 1 / mea.SamplingRate;
			mea.Time =@() (0:length(mea.Data) - 1) * dt - mea.Padding(1);
		end
	case 0
		t = mea.Time();
		assert(abs(1/diff(t(1:2)) - mea.SamplingRate) < tol, 'Time and SamplingRate mismatch.');
end
mea.Duration = t(end) - mea.Padding(2);

end

function mea = check_data(mea, tol)

if ndims(mea.Data) == 3
	[mea.Data, position, map] = reshape_data(mea.Data);

	if isempty(mea.Position)
		mea.Position = position;
	else
		assert(all(abs(mea.Position - position) < tol), 'Check Position.');
	end
	
	if isempty(mea.Map)
		mea.Map = map;
	else
		assert(all(abs(mea.Map - map) < tol), 'Check Map.');
	end
end

mea.Data = int16(rescale(mea.Data, -2^15, 2^15));

if ndims(mea.firing_rate) == 3
    mea.firing_rate = reshape_data(mea.firing_rate);
end
end

function mea = check_position(mea, tol)

switch isempty(mea.Position) + isempty(mea.Map)
	case 0  % Neither is empty - make sure they agree
		map = mea.Map(~isnan(mea.Map));
		[p1, p2] = ind2sub(size(mea.Map), map);
		assert(all(abs([p1; p2] - mea.Position(:)) < tol), 'Position and Map do not match.')
	case 1  % One is empty - use one to make the other
		if isempty(mea.Position)
			map = mea.Map(~isnan(mea.Map));
			[p1, p2] = ind2sub(size(mea.Map), map);
			mea.Position = [p1, p2];
		else
			position = mea.Position;
			mea.Map = nan(max(position));
			inds = sub2ind(max(position), position(:, 1), position(:, 2));
			mea.Map(inds) = 1:length(position);
		end
	case 2  % Both are empty - at least one must be provided for 2D data
		error('''Position'' or ''Map'' must be provided if Data is not 3D.')		
end

end

function [data_mat, position, map] = reshape_data(data_mat)

[nt, time_dim] = max(size(data_mat));

data_mat = shiftdim(data_mat, time_dim - 1);
[~, nr, nc] = size(data_mat);
[xx, yy] = ndgrid(1:nr, 1:nc);
position = [xx(:) yy(:)];
data_mat = reshape(data_mat, nt, nr * nc);

map = reshape(1:nr*nc, nr, nc);
end
