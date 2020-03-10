function [output, mea] = filter_mea(mea, bands, custom_band)
% Filter mea struct to any of the following frequency bands
%     'mua', 'lfp', 'highg'
% Inputs:
%	mea: mea structure containing fields Data, Position, SamplingRate
%	bands: cell or character string indicating which frequency band to
%	       isolate. Default limits for each band are
%				lfp: 2-50 Hz (these data are downsampled to 1 kHz)
%				highg: 50-300 Hz
%               mua: 300-3000 Hz
%	custom_band: 1x2 array indicating custom frequency band. To use a
%		custom band, <band> must be a string (not a cell of strings).
%	Examples: 
%       The following two result in the same filtered data:
%             filter_mea(mea, 'mua');
%             filter_mea(mea, 'mua', [300 3000]);  
%       To filter more than one default band
%			  filter_mea(mea, {'lfp', 'mua'});
%
% Output: Filtered data are saved to <output> struct with fieldname
% matching the <bands>. If <bands='lfp'>, then <skipfactor> for
% downsampling from original data is also saved to <output>

if ~isstruct(mea)
	if ~mea.Properties.Writable
		mea = load(mea.Properties.Source); 
	end
end

if ~isempty(bands) && isnumeric(bands), custom_band = bands; bands = {'custom'}; end

if ~exist('bands', 'var') || isempty(bands); 
	if nargin == 3
		bands = {'custom'};
	else
		bands = {'mua'}; 
	end  % Default to mua if empty
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
end

ElectrodeXY = mea.Position;

if any(strcmpi(bands, 'lfp'))
	% Filter to 2-50 Hz using 150 order bp filter; downsample to 1000 Hz
	disp('Filtering lfp band...')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', SamplingRate);
	skipfactor = max(round(SamplingRate / 1e3), 1);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	temp = downsample(temp, skipfactor);
	output.lfp = temp;
	output.skipfactor = skipfactor;
	if ~isstruct(mea), disp('Writing to file...'), end
	mea.lfp = temp;
	mea.skipfactor = skipfactor;
	disp('Done.')
	clear temp;
end

if any(strcmpi(bands, 'mua'))
	disp('Filtering mua band...')
	b = [3e2, 3e3];
	hipass = min(round(SamplingRate / 2 - 1), b(2));
	if hipass < b(2), warning('Setting CutoffFrequency2 to %f', hipass), end
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',b(1),'CutoffFrequency2',hipass, ...
		'SampleRate',SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	output.mua = temp;
	if ~isstruct(mea), disp('Writing to file...'), end
	mea.mua = temp;
	disp('Done.')
	clear temp;
end

if any(strcmpi(bands, 'custom'))
	disp('Filtering custom band...')
	b = custom_band;
	hipass = min(round(SamplingRate / 2 - 1), b(2));
	if hipass < b(2), warning('Setting CutoffFrequency2 to %f', hipass), end
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',b(1),'CutoffFrequency2',hipass, ...
		'SampleRate',SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	output.custom = temp;
	if ~isstruct(mea), disp('Writing to file...'), end
	mea.custom = temp;
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





