function [delays] = dblr(phif, f, delaytype, Cf, confC)

dim = 1;  % frequency dimension
nanflag = 'includenan';  % don't interpolate nans
method = 'movmean';
df = diff(f(1:2));

switch lower(delaytype)
	case 'phase'
		phif(Cf <= confC) = nan;  % set insignificant values to nan
		delays = smoothdata(phif ./ f', dim, method, winsz / df, nanflag);
	case 'group'
		[~, dphi] = gradient(unwrap(phif), df);  % Compute the gradient of phi wrt freq
		dphi(Cf <= confC) = nan;  % set insignificant values to nan
		delays = smoothdata(dphi, dim, method, winsz / df, nanflag);
	otherwise
		error('delaytype %s unrecognized.', delaytype)
end

end