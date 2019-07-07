% delays_by_localregression
if ~isinteger(phi); phi = single(phi) / 1e4; end
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

for ii = 1:10
    temp = delays(:, :, ii);
    imagesc(t, f, temp, quantile(temp(:), [.05, .95]));
    axis xy; colorbar; ylim([0 100]); drawnow(); pause();
end 

