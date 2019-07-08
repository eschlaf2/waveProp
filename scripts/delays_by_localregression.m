% delays_by_localregression
if isinteger(phi); phi = single(phi) / 1e4; end
phi(C <= confC) = nan;
dphi = permute(phi, [2 1 3]);  % permute array to f x t x pairs
dphi = gradient(dphi, f(2) - f(1));  % d/df

winsz = 3;  % Hz
dim = 1; % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
degree = 2;  % linear
method = 'rlowess';

delays = matlab.internal.math.localRegression(dphi, winsz, dim, ...
            nanflag, degree, method, f);
        
%% Imagesc delays

clims = quantile(delays(:), [.025 .975]);
for ii = 1:10
    temp = delays(:, :, ii);
    imagesc(t, f, temp, clims);
    axis xy; colorbar; ylim([0 500]); drawnow(); pause();
end 

%% Fit delays

predictors = [ones(size(f)); f; f.^2; f.^3]';
delaysR = reshape(delays, length(f), []);
[polyfit, ~, ~, ~, stats] = arrayfun(@(ii) regress(delaysR(:, ii), predictors), 1:length(delaysR), 'uni', 0);
polyfit = reshape(cat(2, polyfit{:}), size(predictors, 2), length(t), []);
polyfit = permute(polyfit, [2 3 1]);
p = reshape(cellfun(@(x) x(3), stats), length(t), []);




