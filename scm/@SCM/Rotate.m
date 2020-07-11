function rotate(P, theta)
    center = P.centerNP;
    rot =@(xy) round( ...
            (xy - center) ...  % translate to center
            * [cos(theta), -sin(theta); sin(theta), cos(theta)] ...  % rotate
        ) + center;  % translate back

    P.stim_center = rot(P.stim_center);

    for ii = 1:size(P.source, 3)
        mat = P.source(:, :, ii);
        rot_mat = false(size(mat));
        [xx, yy] = ind2sub(size(mat), find(mat));
        xnew = rot([xx, yy]);
        inds = sub2ind(P.grid_size, xnew(:, 1), xnew(:, 2));
        rot_mat(inds) = true;
        P.source(:, :, ii) = rot_mat;
    end


end