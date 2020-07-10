function [delays] = dblr(phif, f, delaytype, Cf, confC, winsz)

dim = 1;  % frequency dimension
nanflag = 'omitnan';  % interpolate nans ('includenan' to not interpolate them)
method = 'movmean';
% df = mean(diff(f));
if ~exist('winsz', 'var') || isempty(winsz); winsz = 20; end  % winsz in samples to smooth over
mask = Cf <= confC;

switch lower(delaytype)
	case 'phase'
		phif(mask) = nan;  % set insignificant values to nan
		delays = smoothdata(phif ./ f', dim, method, winsz, ...
			nanflag, 'SamplePoints', f);
	case 'group'
        phif(mask) = nan;  % set insignificant values to nan
% 		[~, dphi] = gradient(unwrap(phif), df);  % Compute the gradient of phi wrt freq
% 		[~, dphi] = phase_gradient(phif);
%		[~, dphi] = gradient(unwrap(phif), df);  % Compute the gradient of phi wrt freq
		dphi = diff_phase(phif);
		df = diff(f);
% 		dphi(Cf <= confC) = nan;
		delays = smoothdata(dphi ./ df(:), dim, method, winsz, ...
			nanflag, 'SamplePoints', f(2:end));
		N = size(delays, 2);
		delays = [nan(1, N); delays];
        delays(mask) = nan;
	otherwise
		error('delaytype %s unrecognized.', delaytype)
end

end
