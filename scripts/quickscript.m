% pat = 'c7'; seizure = 1; T = 10; W = 2; DS = 1e3; units = 1;

if isempty(T), T = 10; else, if ischar(T), T = str2double(T); end, end
if isempty(W), W = 2; else, if ischar(W), W = str2double(W); end, end
if isempty(DS), DS = 1e3; else, if ischar(DS), DS = str2double(DS); end, end
if isempty(units), units = 1; else, if ischar(units), units = str2double(units); end, end  
if ~exist('toi', 'var') || isempty(toi), toi = [-Inf Inf]; end

% units is samples per 1/units sec (i.e. set units=1 for Hz, units=1e3 for
% kHz). Note that if units is 10, for example, then T is in tenths of
% seconds while W is in deca(?)Hz.

basename = compute_coherograms(pat, seizure, T, W, DS, units, toi);

function basename = compute_coherograms(pat, seizure, T, W, DS, units, toi)
%%
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);
time = mea.Time();
if 1.1*T >= diff(toi); toi = (toi - mean(toi)) * 1.1*T + mean(toi); end
mask = time < toi(1) | time > toi(2);
if 2.1*T > diff(toi)
	pad = (2.1*T - diff(toi)) / 2;
	mea.Data(mask, :) = 0;
	mask = time < (toi(1) - pad) | time > (toi(2) + pad);
end
mea.Data(mask, :) = [];
time(mask) = [];

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
basename = strrep(sprintf('%s_cohgram_ds%s_T%02d_W%02d_Hz%d_t%d_%d', ...
	name, ...
	strrep(num2str(DS, '%0.0g'), '+', ''), ...
	round(T), W, units, round(toi)), '-', 'M');
% basename = [name '_cohgram_ds_T' num2str(T, '%02d')];
outfile = matfile(basename, 'writable', true);
% mea = exclude_channels(mea);
% skipfactor = mea.SamplingRate / DS;  % Downsample data to ~DS Hz
[nT, nCh] = size(mea.Data);
[X, Y] = meshgrid(1:nCh, (1:nT) / mea.SamplingRate);
[Xq, Yq] = meshgrid(1:nCh, 1/DS:1/DS:nT/mea.SamplingRate);
data = interp2(X, Y, single(mea.Data), Xq, Yq, 'cubic');

% data = data(1.5*DS:3.5*DS, :);

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

% padding = mea.Padding;  % Store to correct t later
clear mea;  % free up memory

%% Initialize arrays
% ii = length(pairs);
% disp('Initializing arrays with last pair...')
% [C{ii}, phi{ii}, ~, ~, ~, ...
% 		t, f, confC{ii}, phistd{ii}] = ...
% 		cohgramc(data1(:, ii), data2(:, ii), movingwin, params);
% disp('Arrays initialized.')
%%	
% Compute coherence and spectrograms for each pair of channels
% parfor ii = 1:length(pairs) - 1
% 	[C{ii}, phi{ii}, ~, ~, ~, ...
% 		~, ~, confC{ii}, phistd{ii}] = ...
% 		cohgramc(data1(:, ii), data2(:, ii), movingwin, params);
% 	save(sprintf('%s_%d', basename, ii), 'C', 'phi', 
% 	disp(['ii=' num2str(ii)]);
% end

% [C, phi, ~, ~, ~, t, f, confC, ~] = ...
%         cohgramc(...
%             data(:, pairs(1, 1)), data(:, pairs(1, 2)), ...
%             movingwin, params);
        
%%
% Initialize variables
% numpairs = nchoosek(nCh, 2);
numpairs = nCh - 1;
slicesize = 10;
numslices = ceil(numpairs / slicesize);
[C, phi, t, f, confC, S12, S1, S2] = deal(cell(1, numslices));

% p = gcp;
data = smoothdata(single(data));
if any(strcmpi(pat, {'sim', 'waves'}))
	data = data + randn(size(data)) * .01 * diff(quantile(data(:), [.01, .99]));
end

for ii = 1:numslices
    
    disp(ii)  % show progress
    
    i0 = (ii - 1) * slicesize + 1;  % index of starting pair
    iF = min(i0 + slicesize - 1, numpairs);  % index of final pair
    
    tic  % start timer
    [Ct, phit, S12t, S1t, S2t, t{ii}, f{ii}, confC{ii}, ~] = cohgramc(...
            data(:, pairs(i0:iF, 1)), ...  % data1
            data(:, pairs(i0:iF, 2)), ...  % data2
            movingwin, params);  % parameters
    toc  % end timer and display
    
    C{ii} = int16(Ct * 1e4);  % convert to int16
    phi{ii} = int16(phit * 1e4);
    S12{ii} = S12t;
    S1{ii} = S1t(:, :, 1);
    S2{ii} = S2t;
    
    
end
disp('Saving result.')
f = f{1};
t = t{1} / units + time(1);  % correct for padding
confC = int16(confC{1}(1) * 1e4);
C = cat(3, C{:});
phi = cat(3, phi{:});

% Save results
outfile.C = C;
outfile.phi = phi;
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
outfile.S12 = S12;
outfile.S1 = S1;
outfile.S2 = S2;
% outfile.phistd = phistd;

disp('Done.')

plotmean();

% delete(p);

%% Nested plotting function
    function plotmean
        
%         files = dir('*cohgram*ds*.mat');
        
%         close all
%         for file = {files.name}
%             [~, basename, ~] = fileparts(file{:});
%             disp(basename);
%             try
%                 load(basename);
%             catch ME
%                 disp(['Error loading ' basename])
%                 continue
%             end
%             disp('Loaded')
            
            strinfo = strsplit(basename, {'_', '.'});
            pat = strinfo{1};
            seizure = str2double(strinfo{2}(8:end));
            T = str2double(strinfo{end}(2:end));
        
            figure(); fullwidth(true)
            
            mn = single(quantile(C, .9, 3));
            mn(mn <= confC) = nan;
            phimn = single(phi) / 1e4;
            phimn(C <= confC) = nan;
            phistd = std(phimn, [], 3, 'omitnan');
            phimn = mean(phimn, 3, 'omitnan');

            subplot(2,2,1)
            imagesc(t, f, mn'); colorbar; axis xy;
            title(sprintf('%s Seizure %d\nT=%0.1f\nMean coherence', pat, seizure, T))
            xlabel('Time (s)'); ylabel('Frequency (Hz)')
            ylim([0 50])

            subplot(2, 2, 2);
            imagesc(t, f, phimn'); colorbar; axis xy;
            title(sprintf('%s Seizure %d\nT=%0.1f\nMean phi', pat, seizure, T))
            xlabel('Time (s)'); ylabel('Frequency (Hz)')
            ylim([0 50])

    %         print([fid '_mean_phi'], '-dpng')

            subplot(2, 2, 4);
            imagesc(t, f, phistd'); colorbar; axis xy;
            title(sprintf('%s Seizure %d\nT=%0.1f\nstd phi', pat, seizure, T))
            xlabel('Time (s)'); ylabel('Frequency (Hz)')
            ylim([0 50])

            print(basename, '-dpng')

            close(gcf)
%         end
    end

end