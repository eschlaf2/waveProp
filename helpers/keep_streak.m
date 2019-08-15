function [A] = keep_streak(A, minstreak, thresh, dim)
% Keep only the longest streak of consecutive values of A that are above
% threshold (default: nan). Evaluate streaks along dimension dim (default:
% first non-singleton dimension). If streaks are longer than minstreak
% (default: 0), then they are tossed out.

if ~exist('thresh', 'var') || isempty(thresh), thresh = nan; end
if ~exist('dim', 'var') || isempty(dim), dim = find(size(A) > 1, 1); end
if ~exist('minstreak', 'var') || isempty(minstreak), minstreak = 0; end
dim = dim - 1; A = shiftdim(A, dim);
sz = size(A);

A = reshape(A, size(A, 1), []);
cols = size(A, 2);
temp = [nan(1, cols); A; nan(1, cols)];
switches = diff(temp > thresh);
mask = false(size(A));

for ii = 1:cols
    starts = find(switches(:, ii) == 1);
    ends = find(switches(:, ii) == -1) - 1;
    
    [mx, imx] = max(ends - starts);
    if mx < minstreak, continue, end
    on = starts(imx);
    off = ends(imx);
    mask(on:off, ii) = true;
end

A(~mask) = nan;
A = shiftdim(reshape(A, sz), numel(sz) - dim);


