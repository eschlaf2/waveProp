% pat = 'SIM'; seizure = 7; T = 2; W = 2; DS = 1e3; units = 1; type = 'pb';

if isempty(T), T = 10; else, if ischar(T), T = str2double(T); end, end
if isempty(W), W = 2; else, if ischar(W), W = str2double(W); end, end
if isempty(DS), DS = 1e3; else, if ischar(DS), DS = str2double(DS); end, end
if isempty(units), units = 1; else, if ischar(units), units = str2double(units); end, end  
if ~exist('toi', 'var') || isempty(toi), toi = [-Inf Inf]; end
if ~exist('cohfun', 'var') || isempty(cohfun), cohfun = 'c'; end

% units is samples per 1/units sec (i.e. set units=1 for Hz, units=1e3 for
% kHz). Note that if units is 10, for example, then T is in tenths of
% seconds while W is in deca(?)Hz.

basename = compute_coherograms(pat, seizure, T, W, DS, units, toi, cohfun);

function basename = compute_coherograms(pat, seizure, T, W, DS, units, toi, cohfun)
%%
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);
time = mea.Time();
toi(1) = max(toi(1), time(1)); toi(2) = min(toi(2), time(end));

if T >= range(toi); T = floor(range(toi)); fprintf('Using T=%d\n', T); end  % ensure T is not longer than toi
mask = time < toi(1) | time > toi(2);  % Remove datapoints outside of toi
mea.Data(mask, :) = [];
time(mask) = [];

if 2.1*T >= range(time)  % Pad short datasets to ensure multiple time points
    dt = diff(time(1:2));
    pad = (2.1 * T - range(time)) / 2;
    time = time(1) - pad : dt : time(end) + pad;
    if mod(length(time) - length(mea.Data), 2); time = time(1:end-1); end
    mea.Data = padarray(mea.Data, ...  % pad data
        [(length(time) - length(mea.Data)) / 2, 0], ...  % to match length of time
        'both', 'rep');  % by repeating the last value on the beginning and end
end

basename = strrep(sprintf('%s_cohgram_ds%s_T%02d_W%02d_Hz%d_t%d_%d_%s', ...
	name, ...
	strrep(num2str(DS, '%0.0g'), '+', ''), ...
	round(T), W, units, round(toi), cohfun), '-', 'M');
outfile = matfile(basename, 'writable', true);
% mea = exclude_channels(mea);
[nT, nCh] = size(mea.Data);
%%
if strcmpi(cohfun, 'pb')
	
	[~, ~, mea] = mua_events(mea);
	data = false(size(mea.Data));
	data(mea.event_inds) = true;
	binsz = round(mea.SamplingRate / DS);
	nT = floor(nT / binsz);  % new data length
	data = squeeze(sum(reshape(data(1:nT * binsz, :), nT, binsz, nCh), 2));
	time = time(1:binsz:nT*binsz);
	DS = 1 / (diff(time(1:2)));
else
	
	[X, Y] = meshgrid(1:nCh, (1:nT) / mea.SamplingRate);
	[Xq, Yq] = meshgrid(1:nCh, 1/DS:1/DS:nT/mea.SamplingRate);
	data = interp2(X, Y, single(mea.Data), Xq, Yq, 'cubic');
end
data(:, mea.BadChannels) = [];
Fs = DS / units;  % sampling frequency (Hz * units)

nCh = size(data, 2);

%% Set some parameters
STEP = .01 * units;  % Step (s)
THRESH = 5e-3;  % significance threshold

%% Convert parameters to function input
movingwin = [T STEP];  % [window step] seconds
params.err = [1 THRESH];  % [type threshold]
params.Fs = Fs;  % sampling rate (Hz)
params.fpass = [0 100];  % lfp filtered range
params.tapers = [W T 1];  % [bandwidth time k] (numtapers = 2TW - k)
params.pad = max(ceil(log2(20 / T)), 0);  % pad fft filter such that df < .05

position = mea.Position;
badchannels = mea.BadChannels;
position(badchannels, :) = [];  % Remove bad channels
[~, center] = min(sum((position - mean(position)).^2, 2));
pairs = [repmat(center, 1, nCh); 1:nCh]';
pairs(pairs(:, 1) - pairs(:, 2) == 0, :) = [];

% Save information about the run
outfile.data = data;
outfile.position = position;
outfile.badchannels = badchannels;
outfile.units = units;
outfile.pairs = pairs;  % All pairs that include the central electrode
outfile.center = center;
outfile.params = params;
outfile.movingwin = movingwin;
outfile.basename = basename;
outfile.pat = pat;
outfile.seizure = seizure;

% padding = mea.Padding;  % Store to correct t later
clear mea;  % free up memory

%%
% Initialize variables

numpairs = nCh - 1;
slicesize = 10;
numslices = ceil(numpairs / slicesize);
[C, phi, t, f, confC, S12m, S12a, S1, S2] = deal(cell(1, numslices));
if strcmpi(cohfun, 'pb'), zerosp = cell(1, numslices); end
if strcmpi(cohfun, 'c'), data = smoothdata(single(data)); end

if any(strcmpi(pat, {'sim', 'waves'}))
    disp('Adding noise to simulated data')
	if strcmpi(cohfun, 'c')
		data = data + randn(size(data)) * .01 * diff(quantile(data(:), [.01, .99]));
	else
		noise = .2 * randn(size(data));
		data = round(noise) + data; data(data < 0) = 0;
	end
end

%%
for ii = 1:numslices
    
    disp(ii)  % show progress
    
    i0 = (ii - 1) * slicesize + 1;  % index of starting pair
    iF = min(i0 + slicesize - 1, numpairs);  % index of final pair
    
    tic  % start timer
	if strcmpi(cohfun, 'c')
    [Ct, phit, S12t, S1t, S2t, t{ii}, f{ii}, confC{ii}, ~] = cohgramc(...
            data(:, pairs(i0:iF, 1)), ...  % data1
            data(:, pairs(i0:iF, 2)), ...  % data2
            movingwin, params);  % parameters
	else
		[Ct,phit,S12t,S1t,S2t,t{ii},f{ii},zerosp{ii},confC{ii}, ~] = cohgrampb(...
			data(:, pairs(i0:iF, 1)), ...  % data1
            data(:, pairs(i0:iF, 2)), ...  % data2
            movingwin, params);  % parameters
	end
    toc  % end timer and display
    
    
    C{ii} = int16(Ct * 1e4);  % convert to int16
    phi{ii} = int16(phit * 1e4);
    S12m{ii} = int16(1e3 * log10(abs(S12t) ./ max(abs(S12t), 2)));
    S12a{ii} = int16(1e4 * angle(S12t));
    S2{ii} = int16(1e3 * log10(S2t ./ max(S2t, 2)));
    
    
end
disp('Saving result.')
f = f{1};
t = t{1} / units + time(1);  % correct for padding
confC = int16(confC{1}(1) * 1e4);
C = cat(3, C{:});
phi = cat(3, phi{:});
S12m = cat(3, S12m{:});
S12a = cat(3, S12a{:});
S1 = int16(1e3 * log10(S1t(:, :, 1) ./ max(S1t(:, :, 1), 2)));
S2 = cat(3, S2{:});
if strcmpi(cohfun, 'pb'), outfile.zerosp = logical(cat(2, zerosp{:})); end

% Save results
outfile.C = C;
outfile.phi = phi;
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
outfile.S12m = S12m;
outfile.S12a = S12a;
outfile.S1 = S1;
outfile.S2 = S2;
% outfile.phistd = phistd;

disp('Done.')
disp(basename);

end