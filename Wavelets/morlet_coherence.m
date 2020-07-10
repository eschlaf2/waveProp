function [C, phi, S12, S1, S2, freq, scale, coi, confC, phistd] = ...
	morlet_coherence(data1, data2, params)

if nargin < 3 || isempty(params), params = struct; end
[data1, data2] = local_check_data(data1, data2);
P = parse_params(params);
% tapers = dpss(size(data1, 1), P.tapers(1), P.tapers(2)) * sqrt(P.Fs);
args = {P.Fs, P.fpass(1), P.fpass(2), P.omega, 0};
[wave1, period, scale, coi] = local_apply_tapers(data1, tapers, args);
wave2 = local_apply_tapers(data2, tapers, args);

freq = 1 ./ period;

xf =@(data) squeeze(mean(data, 3));
S1 = xf(conj(wave1) .* wave1);
S2 = xf(conj(wave2) .* wave2);
S12 = xf(conj(wave1) .* wave2);
C12 = S12 ./ sqrt(S1 .* S2);
C = abs(C12);
phi = angle(C12);

if nargout > 8
	[confC, phistd] = local_confidence_bounds(wave1, wave2, C, P.err); 
end

end

function  [confC, phistd] = local_confidence_bounds(wave1, wave2, C, err)

for tt = 1:size(wave1, 2)
	J1 = squeeze(wave1(:, tt, :, :));
	J2 = squeeze(wave2(:, tt, :, :));
	N1 = size(J1, 3);
	N2 = size(J2, 3);
	if N1 ~= N2
		if N1 > N2
			J2 = repmat(J2, 1, 1, N1);
		else
			J1 = repmat(J1, 1, 1, N2);
		end
	end
	[confC,phistd]=coherr(squeeze(C(:, tt, :)),J1,J2,err,0);
end

end

function [dataW, period, scale, coi] = ...
	local_apply_tapers(data, tapers, args)

Nd = size(data);
NT = size(tapers);
assert(Nd(1) == NT(1), ...
	'Taper length does not match data length.');

T = tapers(:, :, ones(1, Nd(2)));
D = data(:, :, ones(1, NT(2)));
D = permute(D, [1 3 2]);
data = D.*T;

[dataW, period, scale, coi] = basewave4(data(:, 1, 1), args{:});
dataW = dataW';
dataW(:, :, NT(2), Nd(2)) = nan;

for ii = 1:Nd(2)
	for jj = 1:NT(2)
		dataW(:, :, jj, ii) = basewave4(data(:, jj, ii), args{:})';
	end
end

end

function [data1, data2] = local_check_data(data1, data2)
if all(size(data1) == size(data2)), return, end
data = {data1, data2};
vectors = cellfun(@isvector, data);
assert(any(vectors), 'Sizes of data are mismatched.');
data{vectors} = data{vectors}(:);
assert(size(data{1}, 1) == size(data{2}, 1));
N = cellfun(@(c) size(c, 2), data);
data{vectors} = data{vectors} .* ones(1, max(N)); 
[data1, data2] = data{:}; 


end

function params = parse_params(params)

p = inputParser;
addParameter(p, 'Fs', 1);
addParameter(p, 'fpass', [1 50]);
addParameter(p, 'omega', 6);
addParameter(p, 'err', 0);
addParameter(p, 'tapers', [3 5]);
addParameter(p, 'trialave', 0);
parse(p, params);

params = p.Results;
if numel(params.tapers) == 3
	tapers = params.tapers;
	T = tapers(1);
	W = tapers(2);
	p = tapers(3);
	params.tapers = [T*W 2*T*W-p];
end

end

