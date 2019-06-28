% pat = 'c7'; seizure = 1; T = 10;
compute_coherograms(pat, seizure, T);

function compute_coherograms(pat, seizure, T)

if isempty('T'); T = 10; end
datapath = genpath(['/projectnb/ecog/Data' filesep pat]);  % matlab doesn't follow symlinks so 
addpath(datapath);  % ... add the original data path first
patpath = genpath(pat);  % ... and then add the local patient path on top 
addpath(patpath);  % ... so that it is searched first
% computetimesmethod = 1;
% showplots = false;

fname = sprintf('%s_Seizure%d_Neuroport_10_10.mat', pat, seizure);
mea = load(fname);
[~, name, ~] = fileparts(fname);

% mea = load('SIM/seizing_cortical_field_sim.mat');
% name = mea.Name;
basename = [name '_cohgram_mua_T' num2str(T, '%02d')];
outfile = matfile(basename, 'writable', true);
mea = exclude_channels(mea);

skipfactor = floor(mea.SamplingRate / 1e3);
data = single(downsample(mea.Data, skipfactor));
Fs = mea.SamplingRate / skipfactor;

% [~, mea] = filter_mea(mea, 'mua');
nCh = size(data, 2);
% [~, mea] = get_discharge_times(mea, 'method', computetimesmethod);
% mea.Time = mea.Time();

% T = 10;  % Window (s)
STEP = .5;  % Step (s)
THRESH = 5e-3;  % significance threshold
W = 2;  % bandwidth (Hz)
FS = Fs;  % sampling frequency (Hz)
% FPASS = [0 100];  % Frequencies of interest

movingwin = [T STEP];  % [window step] seconds
params.err = [1 THRESH];  % [type threshold]
params.Fs = FS;  % sampling rate (Hz)
params.fpass = [1 500];  % lfp filtered range
params.tapers = [W T 1]; 

pairs = nchoosek(1:nCh, 2);
data1 = data(:, pairs(:, 1));
data2 = data(:, pairs(:, 2));

outfile.data = data;
clear data;
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

for ii = 0:100:400
    params.fpass = [ii ii+100];
    disp(ii)
    tic
    [C, phi, ~, ~, ~, t, f, confC, ~] = ...
        cohgramc(data1, data2, movingwin, params);
    toc
    outfile.(sprintf('C%03d', ii)) = C;
    outfile.(sprintf('phi%03d', ii)) = phi;
    outfile.(sprintf('f%03d', ii)) = f;
    
    plotmean();  % nested plotting function
    
    clear C phi
end
disp('Saving result.')
t = t - mea.Padding(1);  % correct for padding
confC = confC(1);
% Save results
% outfile.C = C;
% outfile.phi = phi;
outfile.t = t;
% outfile.f = f;
outfile.confC = confC;
% outfile.phistd = phistd;
outfile.position = mea.Position;
outfile.badchannels = mea.BadChannels;
outfile.pairs = pairs;
outfile.params = params;
outfile.movingwin = movingwin;

disp('Done.')

%% Nested plotting function
function plotmean
    disp('Plotting mean')

    mn = mean(C, 3);
    mn(mn < confC(1)) = nan;
    mn(sum(C > confC(1), 3) < size(C, 3) / 2) = nan;  % Set to nan if fewer than half of the values have significant coherence
    figure(); fullwidth(1)
    h = pcolor(t, f, mn'); h.LineStyle = 'none';
    title(sprintf('%s Seizure%d\nT=%f', pat, seizure, T))
    print(gcf, sprintf('%s_%d-%d', basename, round(f(1)), round(f(end))), '-dpng'); 
    close(gcf);
end

end