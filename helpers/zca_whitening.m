function [Y, Z] = zca_whitening(X)
% ZCA whitening of matrix X (M x N)
%	M rows are observations
%	N columns are variables

%%
X = X - mean(X);
sigma = cov(X);  % (N x N)
[U, S, ~] = svd(sigma);  % U * S * V' = sigma

Z = U * 1/sqrt(S + 1e-6) * U';
Y = X * Z;