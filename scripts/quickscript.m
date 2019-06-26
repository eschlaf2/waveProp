% pat = 'c7'; seizure = 1; T = 10;
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first
computetimesmethod = 1;
showplots = false;

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
basename = [name '_cohgram_T' num2str(T, '%02d')];
outfile = matfile(basename, 'writable', true);
% mea = exclude_channels(mea);
mea = exclude_channels(mea);
[~, mea] = filter_mea(mea, 'lfp');
nCh = size(mea.lfp, 2);
% [~, mea] = get_discharge_times(mea, 'method', computetimesmethod);
% mea.Time = mea.Time();

% T = 10;  % Window (s)
STEP = .5;  % Step (s)
THRESH = 5e-2;  % significance threshold
TW = 20;  % bandwidth (Hz)
FS = floor(mea.SamplingRate / mea.skipfactor);  % sampling frequency (Hz)
% FPASS = [0 100];  % Frequencies of interest

movingwin = [T STEP];  % [window step] seconds
params.err = [1 THRESH];  % [type threshold]
params.Fs = FS;  % sampling rate (Hz)
params.fpass = [2 50];  % lfp filtered range
params.tapers = [TW 2*TW-1]; 

pairs = nchoosek(1:nCh, 2);
data1 = mea.lfp(:, pairs(:, 1));
data2 = mea.lfp(:, pairs(:, 2));

%% Initialize arrays
ii = length(pairs);
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

[C, phi, ~, ~, ~, t, f, confC, phistd] = ...
	cohgramc(data1, data2, movingwin, params);

disp('Saving result.')
t = t - mea.Padding(1);  % correct for padding
% Save results
outfile.C = C;
outfile.phi = phi;
outfile.t = t;
outfile.f = f;
outfile.confC = confC;
outfile.phistd = phistd;
outfile.position = mea.Position;
outfile.badchannels = mea.BadChannels;
outfile.pairs = pairs;
outfile.params = params;
outfile.lfp = mea.lfp;
outfile.movingwin = movingwin;

disp('Done.')

