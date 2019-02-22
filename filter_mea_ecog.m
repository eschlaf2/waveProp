function [mea, ecog] = filter_mea_ecog(mea, ecog, outfile, bands)


MEA = exist('mea', 'var') && ~isempty(mea);
ECOG = exist('ecog', 'var') && ~isempty(ecog);
if ~exist('outfile', 'var') || isempty(outfile)
	if MEA 
		outfile = [mea.Name '_Filt'];
	end
	if ECOG
		outfile = [ecog.Name '_Filt'];
	end
end
if ~exist('bands', 'var') || isempty(bands)
	bands = {'mua'};
end


if MEA
	
	if any(strcmpi(fieldnames(mea), 'Position'))  % deprecated struct field
		mea.ElectrodeXY = mea.Position;
		% mea.Position = 'Renamed as ElectrodeXY';
	end

	data = mea.Data;
	SamplingRate = mea.SamplingRate;
	BadChannels = mea.BadChannels;
	
	ElectrodeXY = mea.ElectrodeXY;

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
end

if ECOG

	data = ecog.Data;
	SamplingRate = ecog.SamplingRate;
	BadChannels = ecog.BadChannels;
	ElectrodeXY = ecog.ElectrodeXY;

	if isstruct(ecog)
		
		if isfield(ecog, 'Position')  % deprecated
			ecog.ElectrodeXY = ecog.Position;
			rmfield(ecog, 'Position')
		end

		rmfield(ecog, 'Data');
		disp('Converting ecog to matfile.')
		save(outfile, '-v7.3', '-struct', 'ecog');
		clear ecog
		ecog = matfile(outfile, 'writable', true)
	end

	disp('Filtering ecog to lfp band.')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	ecog.lfp = temp

	disp('Filtering ecog to high-gamma band.')
	bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
		'CutoffFrequency1', 50, 'CutoffFrequency2', 124, ...  % 150 in the Schevon paper, but max is 124.93 given sampling freq
		'SampleRate', SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	ecog.highg = temp;
	
	ecog.X = ElectrodeXY(:, 1);
	ecog.X(BadChannels) = [];
	ecog.Y = ElectrodeXY(:, 2);
	ecog.Y(BadChannels) = [];
end




