function sources = count_wave_sources(wave_fits)

N_edges = 65;
edges = linspace(-pi, pi, N_edges);
halfwin = 5;  % (s)
tt = wave_fits.computeTimes / 1e3;
time = tt(1):.1:tt(end);
% Rotate Z so that the distributions are centered and it is easy to see
% deviations
% 			rotateby = angle(sum(exp(1j*res(iF).Z), 'omitnan'));
rotateby = 0;
Z = angle(exp(1j * (wave_fits.Z - rotateby)));
N = nan(numel(time), N_edges - 1);

for tidx = numel(time):-1:1
	mask = (tt >= (time(tidx) - halfwin)) & (tt <= (time(tidx) + halfwin));
	N(tidx, :) = histcounts(Z(mask), edges, 'normalization', 'probability');
end

% win = gausswin(round(2*halfwin / diff(time(1:2)))) * ...
% 	gausswin(round(pi / 2 / diff(edges(1:2))))';  % Blur with gaussian kernel: 2*halfwin by pi/2
win = gausswin(round(pi / 2 / diff(edges(1:2))))';
Zdens = conv2(repmat(N(:, :), 1, 3), win, 'same');  % Circular padding and convolution
Zdens = Zdens(:, N_edges-1:2*(N_edges-1));  % Take the middle of the result
MZ = max(Zdens);
MZ = MZ / sum(MZ);
[pks, locs, w, p] = findpeaks(MZ);

% peaks_mask = diff(circshift(MZ, 1)) > 0 & diff(MZ) <= 0 & zscore(MZ(1:end-1)) > 2;

sources.max_density = MZ;
sources.edges = edges;
sources.peaks = pks;
sources.locs = locs;
sources.width = w;
sources.prominence = p;


