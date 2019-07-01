% pat = 'c7'; seizure = 1; T = 10;
compute_coherograms(pat, seizure, T);

function compute_coherograms(pat, seizure, T)

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
basename = [name '_cohgram_ds_T' num2str(T, '%02d')];
outfile = matfile(basename, 'writable', true);
mea = exclude_channels(mea);

skipfactor = floor(mea.SamplingRate / 1e3);  % Downsample data to ~1e3 Hz
data = downsample(mea.Data, skipfactor);
Fs = mea.SamplingRate / skipfactor;

nCh = size(data, 2);

%% Set some parameters
% T = 10;  % Window (s)
STEP = .5;  % Step (s)
THRESH = 5e-3;  % significance threshold
W = 2;  % bandwidth (Hz)
FS = Fs;  % sampling frequency (Hz)
% FPASS = [0 100];  % Frequencies of interest

%% Convert parameters to function input
movingwin = [T STEP];  % [window step] seconds
params.err = [1 THRESH];  % [type threshold]
params.Fs = FS;  % sampling rate (Hz)
params.fpass = [0 500];  % lfp filtered range
params.tapers = [W T 1];  % [bandwidth time k] (numtapers = 2TW - k)
params.pad = -1;  % no padding

mea.Position(mea.BadChannels, :) = [];  % Remove bad channels

% Save information about the run
outfile.data = data;
outfile.position = mea.Position;
outfile.badchannels = mea.BadChannels;
outfile.pairs = nchoosek(1:nCh, 2);  % generate all pairs of channels
outfile.params = params;
outfile.movingwin = movingwin;

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
numpairs = nchoosek(nCh, 2);
slicesize = 100;
numslices = floor(numpairs / slicesize);
[C, phi, t, f, confC] = deal(cell(1, numslices));

parfor ii = 1:numslices
    
    disp(ii)  % show progress
    
    pairs = nchoosek(1:nCh, 2);  % avoid overhead communication
    i0 = (ii - 1) * slicesize + 1;  % index of starting pair
    iF = min(i0 + slicesize - 1, numpairs);  % index of final pair
    
    tic  % start timer
    [C{ii}, phi{ii}, ~, ~, ~, t{ii}, f{ii}, confC{ii}, ~] = cohgramc(...
            single(data(:, pairs(i0:iF, 1))), ...  % data1
            single(data(:, pairs(i0:iF, 2))), ...  % data2
            movingwin, params);  % parameters
    toc  % end timer and display
    
    C{ii} = int16(C{ii} * 1e4);  % convert to int16
    phi{ii} = int16(phi{ii} * 1e4);
    
%     if ii == 1
%         C = int16(Csub * 1e4);
%         phi = int16(phisub * 1e4);
%     else 
%         C = cat(3, C, int16(Csub * 1e4));  % scale C up and cast to int16
%         phi = cat(3, phi, int16(phisub * 1e4));  % scale phi up and then cut it to int16
%     end
%     outfile.(sprintf('C%03d', ii)) = int16(C * 1e4);
%     outfile.(sprintf('phi%03d', ii)) = int16(phi / pi * 1e4);
%     outfile.(sprintf('f%03d', ii)) = f;
        
%     clear Csub phisub
    
end
disp('Saving result.')
plotmean();  % nested plotting function
f = f{1};
t = t{1} - padding(1);  % correct for padding
confC = int16(confC{1}(1) * 1e4);
% Save results
outfile.C = cat(3, C{:});
outfile.phi = cat(3, phi{:});
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
% outfile.phistd = phistd;

disp('Done.')

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
%         C = single(C);
        mn = quantile(C, .9, 3);
%         C(C <= confC) = nan;
        mn(mn < cast(confC, 'like', mn)) = nan;
%         mn(sum(C > confC, 3) < length(pairs) / 2) = nan;
        
        imagesc(t, f, mn'); colorbar; axis xy;
        title(sprintf('%s Seizure %d\nT=%0.1f', pat, seizure, T))
        xlabel('Time (s)'); ylabel('Frequency (Hz)')
        
        print(fid, '-dpng')
        close(gcf)
%         end
    end

end