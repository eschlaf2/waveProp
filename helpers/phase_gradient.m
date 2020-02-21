function [fx, fy] = phase_gradient(X)

X = padarray(exp(1j*X), [1 1], 'replicate');
A = X(:, 1:end-2);
B = X(:, 2:end - 1);
C = X(:, 3:end);
fx = angle(A.*conj(B) + B.*conj(C));
fx = fx(2:end-1, :);

A = X(1:end-2, :);
B = X(2:end-1, :);
C = X(3:end, :);
fy = angle(A.*conj(B) + B.*conj(C));
fy = fy(:, 2:end-1);