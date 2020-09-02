function mask = make_mask(grid_size, shape, boundary_type)
    
    if nargin < 2 || isempty(shape), shape = 'bean'; end
    if nargin < 3 || isempty(boundary_type), boundary_type = 1; end
    shape = validatestring(shape, {'bean', 'disk', 'circle'});
        
    switch shape
        case 'bean'
            [xx, yy] = ndgrid((1:grid_size(1)) - 25 - 15, (1:grid_size(2)) - 25);
            theta = angle(xx + 1j.*yy) - 3*pi/4;
            mask = sqrt(xx.^2 + yy.^2) < 32*(sin(theta).^3 + cos(theta).^3);
            
        case {'disk', 'circle'}
            
            [xx, yy] = ndgrid(1:grid_size(1), 1:grid_size(2));
            R = grid_size ./ (grid_size - 6);
            xx = rescale(xx, -R(1), R(1));
            yy = rescale(yy, -R(2), R(2));
            mask = xx.^2 + yy.^2 <= 1;
            
    end
    
    mask = conv2(mask, ones(3)/9, 'same');
    
    switch boundary_type
        case 1  % mean along boundary
            
        case 2  % 1 / mean along boundary
            mask = 1./mask;
            mask(isinf(mask)) = 0;
        case 3  % strong boundary
            mask(mask > 0 & mask < 1) = 10;
    end
    
    
end