function [ coh, phi, freq, coh_conf ] = compute_coherence( data, params, varargin )
%COMPUTE_COHERENCE Computes the coherence between electrodes using the
%multitaper method.
%   [COH,PHI,FREQ,COH_CONF]=COMPUTE_COHERENCE(DATA,PARAMS) returns the
%   coherence COH (NxNxF) and its phase PHI (NxNxF) for each pair of the N
%   electrodes available in DATA (TxN) at each frequency FREQ (F) where T
%   is the number of time samples. The confidence level at which the
%   coherence is signficantly higher than 0 is also returned in COH_CONF.
%   Parameters used for this computation are all defined in PARAMS (see the
%   help for function coherencyc in Chronux).

nb = size(data,1);
nfft = max(2^(nextpow2(nb) + params.pad), nb);
freq = getfgrid(params.Fs, nfft, params.fpass); 


coh = nan([size(data,2) size(data,2) length(freq)]);
phi = nan([size(data,2) size(data,2) length(freq)]);
coh_conf = nan([size(data,2) size(data,2) length(freq)]);


[X, Y] = ndgrid(1:size(data, 2), 1:size(data, 2));
X = triu(X); X(X == 0) = [];
Y = triu(Y); Y(Y == 0) = [];
pairs = [X(:) Y(:)];
pairs(X(:) == Y(:), :) = [];

for i = 1:numel(varargin)
	switch varargin{i}
		case 'pairs'
			[X, Y] = ndgrid(varargin{i+1}, 1:size(data, 2));
			
% 			coh = nan([size(X) length(freq)]);
% 			phi = nan([size(data,2) size(data,2) length(freq)]);
% 			coh_conf = nan([size(data,2) size(data,2) length(freq)]);
			
			X = triu(X); X(X == 0) = [];
			Y = triu(Y); Y(Y == 0) = [];
			pairs = [X(:) Y(:)];
			pairs(X(:) == Y(:), :) = [];
			
		otherwise
			if ischar(varargin{i})
				warning('Input not recognized: %s', varargin{i});
			end
	end
end

for p = 1:length(pairs)
% for i = 1 : size(data,2)
	i = pairs(p, 1); j = pairs(p, 2);
    d1 = data(:,i);
    d1 = d1 - mean(d1);
% 	percent_prog(p, 1, size(pairs, 1), 'Coherence: ')
%     for j = i+1 : size(data,2)
%         fprintf('Coherence electrodes %d - %d...\n', i, j);
	d2 = data(:,j);
	d2 = d2 - mean(d2);
	[coh(i,j,:), phi(i,j,:), ~, ~, ~, ~, coh_conf(i,j,:), ~] = coherencyc(d1, d2, params);
	coh(j,i,:) = coh(i,j,:);
	coh_conf(j,i,:) = coh_conf(i,j,:);
	phi(j,i,:) = -phi(i,j,:);
%     end
end

end

