function map = DefaultExcitabilityMap(P)
    [ix, iy] = ind2sub(P.grid_size, 1:prod(P.grid_size));
    R = min(P.grid_size) / 2 - 3;
    map = zeros(P.grid_size);
    map(sum(([ix' iy'] - P.grid_size/2).^2, 2) < R.^2) = .5;
end