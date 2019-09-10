function Z = smoothZ(Z)
% smoothZ = @(Z) mod(smoothdata(unwrap(2*Z, 1.2*pi)/2) + pi, 2*pi) - pi;
Z = mod(smoothdata(unwrap(2*Z, 1.2*pi)/2) + pi, 2*pi) - pi;