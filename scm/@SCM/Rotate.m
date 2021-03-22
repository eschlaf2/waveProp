function Rotate(scm, theta)
    center = scm.grid_size/2;
    rot =@(xy) round( ...
            (xy - center) ...  % translate to center
            * [cos(theta), -sin(theta); sin(theta), cos(theta)] ...  % rotate
        ) + center;  % translate back

    scm.stim_center = rot(scm.stim_center);

    for ii = 1:size(scm.source, 3)
        mat = scm.source(:, :, ii);
        
        [Y, X] = find(true(size(mat)));
        xnew = rot([X, Y]);
        F = scatteredInterpolant(X, Y, mat(:), 'linear', 'nearest');
        rot_mat = reshape(F(xnew(:, 1), xnew(:, 2)), scm.grid_size);
        scm.source(:, :, ii) = rot_mat;
    end


end