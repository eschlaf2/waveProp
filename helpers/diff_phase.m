function dphi = diff_phase(phi)

if isreal(phi), phi = exp(1j*phi); end
phi = [phi(1, :); phi; phi(end, :)];
A = phi(1:end-2, :);
B = phi(2:end-1, :);
C = phi(3:end, :);
dphi = angle(A.*conj(B) + B.*conj(C));
dphi = dphi(2:end-1, :);
end