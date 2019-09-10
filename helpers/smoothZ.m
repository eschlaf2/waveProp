function Z = smoothZ(Z)

Z = mod(smoothdata(unwrap(2*Z)/2) + pi, 2*pi) - pi;