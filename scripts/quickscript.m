% pat = 'c7'; seizure = 1; T = 10; W = 2; DS = 1e3;

if isempty(T), T = 10; else, if ischar(T), T = str2double(T); end, end
if isempty(W), W = 2; else, if ischar(W), W = str2double(W); end, end
if isempty(DS), DS = 1e3; else, if ischar(DS), DS = str2double(DS); end, end
compute_coherograms(pat, seizure, T, W);

function compute_coherograms(pat, seizure, T, W, DS)

if isempty(T); T = 10; end
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
basename = sprintf('%s_cohgram_ds_T%02d_W%02d', name, T, W);
% basename = [name '_cohgram_ds_T' num2str(T, '%02d')];
outfile = matfile(basename, 'writable', true);
% mea = exclude_channels(mea);

skipfactor = floor(mea.SamplingRate / DS);  % Downsample data to ~DS Hz
data = downsample(mea.Data, skipfactor);
data(:, mea.BadChannels) = [];
Fs = mea.SamplingRate / skipfactor;

nCh = size(data, 2);

%% Set some parameters
% T = 10;  % Window (s)
STEP = .05;  % Step (s)
THRESH = 5e-4;  % significance threshold
% W = 2;  % bandwidth (Hz)
FS = Fs;  % sampling frequency (Hz)
% FPASS = [0 100];  % Frequencies of interest

%% Convert parameters to function input
movingwin = [T STEP];  % [window step] seconds
params.err = [1 THRESH];  % [type threshold]
params.Fs = FS;  % sampling rate (Hz)
params.fpass = [0 300];  % lfp filtered range
params.tapers = [W T 1];  % [bandwidth time k] (numtapers = 2TW - k)
% params.pad = -1;  % no padding

mea.Position(mea.BadChannels, :) = [];  % Remove bad channels
[~, center] = min(sum((mea.Position - mean(mea.Position)).^2, 2));
pairs = [repmat(center, 1, nCh); 1:nCh]';
pairs(pairs(:, 1) - pairs(:, 2) == 0, :) = [];

% Save information about the run
outfile.data = data;
outfile.position = mea.Position;
outfile.badchannels = mea.BadChannels;
% outfile.pairs = nchoosek(1:nCh, 2);  % generate all pairs of channels
outfile.pairs = pairs;  % All pairs that include the central electrode
outfile.center = center;
outfile.params = params;
outfile.movingwin = movingwin;
outfile.basename = basename;

padding = mea.Padding;  % Store to correct t later
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
[C, phi, t, f, confC] = deal(cell(1, numslices));

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
    [Ct, phit, ~, ~, ~, t{ii}, f{ii}, confC{ii}, ~] = cohgramc(...
            data(:, pairs(i0:iF, 1)), ...  % data1
            data(:, pairs(i0:iF, 2)), ...  % data2
            movingwin, params);  % parameters
    toc  % end timer and display
    
    C{ii} = int16(Ct * 1e4);  % convert to int16
    phi{ii} = int16(phit * 1e4);
    
    
end
disp('Saving result.')
f = f{1};
t = t{1} - padding(1);  % correct for padding
confC = int16(confC{1}(1) * 1e4);
C = cat(3, C{:});
phi = cat(3, phi{:});

% Save results
outfile.C = C;
outfile.phi = phi;
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
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
            
            strinfo = strsplit(basename, '_');
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