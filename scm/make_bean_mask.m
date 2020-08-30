function mask = make_bean_mask(grid_size)
        
    [xx, yy] = ndgrid((1:grid_size(1)) - 25 - 15, (1:grid_size(2)) - 25);
    theta = angle(xx + 1j.*yy) - 3*pi/4;
    mask = sqrt(xx.^2 + yy.^2) < 32*(sin(theta).^3 + cos(theta).^3);
    mask = conv2(mask, ones(3)/9, 'same');

end