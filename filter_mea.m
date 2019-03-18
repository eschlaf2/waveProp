function [mea] = filter_mea(mea, outfile, bands)


if ~exist('outfile', 'var') || isempty(outfile)
	outfile = [mea.Name '_Filt'];
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
if isprop(mea, 'BadChannels')
	BadChannels = mea.BadChannels;
elseif isprop(mea, 'Exclude')
	BadChannels = mea.Exclude;
end

ElectrodeXY = mea.Position;

if isstruct(mea)

	disp('Converting mea to matfile...')
	if ~exist(outfile, 'file')
		save(outfile, '-v7.3', '-struct', 'mea');
	end
	clear mea
	mea = matfile(outfile, 'writable', true);
	disp('Done.')
end

if any(strcmpi(bands, 'lfp'))
	disp('Filtering lfp band...')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', SamplingRate);

	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	disp('Writing to file...')
	mea.lfp = temp;
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
	disp('Writing to file...')
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
	disp('Writing to file...')
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





