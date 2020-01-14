function [coeff, score, explained] = pca_on_windows(windows)
% Computes PCA on different windows of time from a trace. I.e. assumes that
% variables are time relative to the window and observations are different
% windows of time from an original long recording.
% Input:
%	windows: (t x w array) where t is the relative time dimension (find
%		coefficients along this dimension)
% Returns:
%   coeff: (t x t array) PC coefficients
%   score: (w x t array) windows in PC space
%   explained: (t x 1 vector) variance explained by each PC

[coeff, score, ~, ~, explained] = pca(windows');