function [] = create_epoch(pat, seizures, varargin)
%% Function to create an Epoch file

%% Parse input and set defaults
p = inputParser;

addRequired(p, 'pat', @ischar);
addOptional(p, 'seizures', 1:100, @isnumeric);
addParameter(p, 'padding', [60 60], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'datapath', pwd, @ischar);

parse(p, pat, seizures, varargin{:});
struct2var(p.Results)


%% Get metadata and seizures
patPath = genpath(fullfile(datapath, pat));
addpath(patPath);

metadata = jsondecode(fileread([pat '_metadata.json']));
num_seizures = numel(metadata.seizures);

if ~exist('seizures', 'var') || isempty(seizures)
	seizures = 1:num_seizures;
else
	seizures(seizures > num_seizures) = [];
	seizures(seizures < 1) = [];
end

%% Extract each seizure
for s = seizures  % which seizures
    try 
        seizure = metadata.seizures{s}; 
    catch ME
        seizure = metadata.seizures(s); 
    end
    rawFile = seizure.ictal_micros;
    
    nsx = openNSx(rawFile, 'read');
    nsx.Name = sprintf('%s_Seizure%d', pat, seizure.number);
    disp(nsx.Name);
    
    nsx.RawFile = seizure.ictal_micros;
	outMat = sprintf('%s_Seizure%d_Neuroport_%d_%d.mat', ...
		fullfile(datapath, pat, pat), seizure.number, padding);
    nsx.Map = metadata.electrodemap;
    nsx.StartTime = seizure.onset;
    nsx.EndTime = seizure.offset;
    nsx.Padding = padding;  
    nsx.SamplingRate = nsx.MetaTags.SamplingFreq;
	if isfield(seizure, 'exclude')
		nsx.BadChannels = seizure.exclude;
	else
		nsx.BadChannels = [];
	end
%   m = matfile(outMat, 'writable', true);
    
    nsx.PauseInds = nsx.MetaTags.DataPoints;

    % Get position
    [ids, sortOrder] = sort(nsx.Map(:));    
    mask = arrayfun(@(ch) any(ch == ids), nsx.MetaTags.ChannelID);
    nsx.NChannels = nsx.MetaTags.ChannelCount;
    nsx.Position = nan(nsx.NChannels, 2);
    [X, Y] = ndgrid(1:size(nsx.Map, 1), 1:size(nsx.Map, 2));
    sortOrder(ids == -1) = []; ids(ids == -1) = [];
    XY = [X(:) Y(:)];
    XY = XY(sortOrder, :);
    nsx.Position = XY(nsx.MetaTags.ChannelID(mask), :);
    nsx.Labels = nsx.MetaTags.ChannelID(mask);
        
    % Get data
    if nsx.RawData.PausedFile
        nsx.Data = cell2mat(nsx.Data(logical(nsx.MetaTags.Timestamp)));
    end 
    
    [nCh, nSamp] = size(nsx.Data);
    
    if nsx.StartTime < nsx.Padding(1)
        startInd = 1; nsx.Padding(1) = nsx.StartTime;
    else 
        startInd = round((nsx.StartTime - nsx.Padding(1)) * nsx.SamplingRate);
    end
    endInd = (nsx.EndTime + nsx.Padding(2)) * nsx.SamplingRate;
    if endInd > nSamp
        endInd = nSamp;		
        nsx.Padding(2) = (nSamp - (nsx.SamplingRate * nsx.EndTime)) / nsx.SamplingRate;
    end

    nsx.Data = nsx.Data(mask, startInd:endInd);
    nsx.Data = nsx.Data';
    nsx.Duration = nsx.EndTime - nsx.StartTime;
    nsx.Time = @() linspace(-nsx.Padding(1), nsx.Duration + nsx.Padding(2), size(nsx.Data, 1));
	try
		nsx = rmfield(nsx, {'ElectrodesInfo', 'MetaTags', 'RawData'});
	catch ME
		nsx = rmfield(nsx, {'MetaTags', 'RawData'});
	end

    save(outMat, '-v7.3', '-struct', 'nsx');
    clear nsx
    
end
