function [output, mea] = filter_mea(mea, ~, bands)

if ~isstruct(mea)
	if ~mea.Properties.Writable
		mea = load(mea.Properties.Source); 
	end
end

if ~exist('bands', 'var') || isempty(bands)
	bands = {'mua'};
end


if any(strcmpi(fieldnames(mea), 'ElectrodeXY'))  % deprecated struct field
	mea.Position = mea.ElectrodeXY;
	% mea.ElectrodeXY = 'Renamed as Position';
end

data = mea.Data;
SamplingRate = mea.SamplingRate;
BadChannels = [];
if any(strcmpi(fieldnames(mea), 'BadChannels'))
	BadChannels = mea.BadChannels;
elseif any(strcmpi(fieldnames(mea), 'Exclude'))
	BadChannels = mea.Exclude;
	mea.BadChannels = BadChannels;
end

ElectrodeXY = mea.Position;

if any(strcmpi(bands, 'lfp'))
	% Filter to 2-50 Hz using 150 order bp filter; downsample to 1000 Hz
	disp('Filtering lfp band...')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', SamplingRate);
	skipfactor = round(SamplingRate / 1e3);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	temp = downsample(temp, skipfactor);
	output.lfp = temp;
	output.skipfactor = skipfactor;
	if ~isstruct(mea), disp('Writing to file...'), end
	disp('Done.')
	clear temp;
end

if any(strcmpi(bands, 'mua'))
	disp('Filtering mua band...')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',3e2,'CutoffFrequency2',3e3, ...
		'SampleRate',SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	output.mua = temp;
	if ~isstruct(mea), disp('Writing to file...'), end
	mea.mua = temp;
	disp('Done.')
	clear temp;
end

if any(strcmpi(bands, 'highg'))
	disp('Filtering high-gamma band...')
	bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
		'CutoffFrequency1', 50, 'CutoffFrequency2', 300, ...
		'SampleRate', SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	output.highg = temp;
	if ~isstruct(mea), disp('Writing to file...'); end
	mea.highg = temp;
	disp('Done.')
	clear temp;
end

X = ElectrodeXY(:, 1);
X(BadChannels) = [];
mea.X = X;

Y = ElectrodeXY(:, 2);
Y(BadChannels) = [];
mea.Y = Y;





