function [mea, ecog] = filter_mea_ecog(mea, ecog)


MEA = exist('mea', 'var');
ECOG = exist('ecog', 'var');

if MEA

	if isfield(mea, 'Position')  % deprecated
		mea.ElectrodeXY = mea.Position;
		rmfield(mea, 'Position')
	end

	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', mea.SamplingRate);

	mea.lfp = single(filtfilt(bpFilt, double(mea.Data)));
	mea.lfp(:, mea.BadChannels) = [];
	
	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',3e2,'CutoffFrequency2',3e3, ...
		'SampleRate',mea.SamplingRate);
	mea.mua = single(filtfilt(bpFilt, double(mea.Data)));
	mea.mua(:, mea.BadChannels) = [];

	bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
		'CutoffFrequency1', 50, 'CutoffFrequency2', 300, ...
		'SampleRate', mea.SamplingRate);
	mea.highg = single(filtfilt(bpFilt, double(mea.Data)));
	mea.highg(:, mea.BadChannels) = [];
	
	mea.X = mea.ElectrodeXY(:, 1);
	mea.X(mea.BadChannels) = [];

	mea.Y = mea.ElectrodeXY(:, 2);
	mea.Y(mea.BadChannels) = [];
end

if ECOG

	if isfield(ecog, 'Position')  % deprecated
		ecog.ElectrodeXY = ecog.Position;
		rmfield(ecog, 'Position')
	end

	bpFilt = designfilt('bandpassfir','FilterOrder',150, ...
		'CutoffFrequency1',2,'CutoffFrequency2',50, ...
		'SampleRate', ecog.SamplingRate);
	ecog.lfp = single(filtfilt(bpFilt, double(ecog.Data)));
	ecog.lfp(:, ecog.BadChannels) = [];

	bpFilt = designfilt('bandpassfir', 'FilterOrder', 150, ...
		'CutoffFrequency1', 50, 'CutoffFrequency2', 124, ...  % 150 in the Schevon paper, but max is 124.93 given sampling freq
		'SampleRate', ecog.SamplingRate);
	ecog.highg = single(filtfilt(bpFilt, double(ecog.Data)));
	ecog.highg(:, ecog.BadChannels) = [];
	
	ecog.X = ecog.ElectrodeXY(:, 1);
	ecog.X(ecog.BadChannels) = [];
	ecog.Y = ecog.ElectrodeXY(:, 2);
	ecog.Y(ecog.BadChannels) = [];
end




