function [mi, mi_mat] = moran_index(A, pos, map)
    % Moran's index (Martinet, L.-E., 2015). A measure of spatial
    % autocorrelation. If <map> is not given, assume adjacent indicator
    % map. 
    
    % Convert A into a vector and remove nan values
    assert(numel(A) == size(pos, 1), ...
        '<A> and <pos> have different sizes.');
    A = A(:);
    pos(isnan(A), :) = [];
    A(isnan(A)) = [];
    A = A - mean(A);
    N = numel(A);
    T = A * A';
    
    % Default map is adjacent indicator map
    if nargin < 3
        xx = pos(:, 1); yy = pos(:, 2);
        map = zeros(size(T));
        map( (abs(xx - xx') == 1) & (yy == yy') ) = 1;
        map( (abs(yy - yy') == 1) & (xx == xx') ) = 1;
        map = map ./ sum(map, 2);
        map(isnan(map)) = 0;
    end
    assert(all(size(map) == [N, N]));
    
    % Compute mi
    mi = N / sum(map, 'all') * sum(T .* map, 'all') ./ sum(diag(T));
    if nargout > 1, mi_mat = sum(T .* map, 2) ./ sum(diag(T)); end

end