function map = expanding_ring(scm, t)
    
    WIDTH = 2 / scm.dx;  % ring width (in scm spatial units)
    
    r1 = scm.expansion_rate * t;
    r2 = r1 - WIDTH;
    
    gs = scm.grid_size;
    sc = scm.stim_center;
    [xx, yy] = ndgrid((1:gs(1)) - sc(1), (1:gs(2)) - sc(2));
    dist = sqrt(xx.^2 + yy.^2);
    map = false(gs);
    map(r2 <= dist & r1 >= dist) = true;