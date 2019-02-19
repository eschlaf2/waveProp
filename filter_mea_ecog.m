function [mea, ecog] = filter_mea_ecog(mea, ecog, outfile)


MEA = exist('mea', 'var');
ECOG = exist('ecog', 'var');
if ~exist('outfile', 'var')
	if MEA 
		outfile = [mea.Name '_meaFilt'];
	end
	if ECOG
		outfile = [ecog.Name '_ecogFilt'];
	end
end


if MEA
	data = mea.Data;
	SamplingRate = mea.SamplingRate;
	BadChannels = mea.BadChannels;
	ElectrodeXY = mea.ElectrodeXY;

	if isstruct(mea)

		if isfield(mea, 'Position')  % deprecated
			mea.ElectrodeXY = mea.Position;
			rmfield(mea, 'Position');
		end

		rmfield(mea, 'Data');
		disp('Converting mea to matfile.')
		save(outfile, '-v7.3', '-struct', 'mea');
		clear mea
		mea = matfile(outfile, 'writable', true)
	end

	disp('Filtering lfp band.')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', SamplingRate);

	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	mea.lfp = temp;
	clear temp;
	
	disp('Filtering mua band.')
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',3e2,'CutoffFrequency2',3e3, ...
		'SampleRate',SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	mea.mua = temp;
	clear temp;

	disp('Filtering high-gamma band.')
	bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
		'CutoffFrequency1', 50, 'CutoffFrequency2', 300, ...
		'SampleRate', SamplingRate);
	temp = single(filtfilt(bpFilt, double(data)));
	temp(:, BadChannels) = [];
	mea.highg = temp;
	clear temp;
	
	mea.X = ElectrodeXY(:, 1);
	mea.X(BadChannels) = [];

	mea.Y = ElectrodeXY(:, 2);
	mea.Y(BadChannels) = [];
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




