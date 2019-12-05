function output = test_for_recruitment(mea, method, PLOT)
% Test to see if the electrode array was recruited to the seizure as
% described in ?Schevon, Catherine A., et al. ?Evidence of an Inhibitory
% Restraint of Seizure Activity in Humans.? Nature Communications, vol. 3,
% no. 1, Nature Publishing Group, Jan. 2012, p. 1060,
% doi:10.1038/ncomms2056.
% 
% Uses weighted phase locking index

METHOD_ALL = {'fano', 'wpli'};
if ~exist('method', 'var') || isempty(method)
	method = {'fano', 'wpli'};
else
	switch class(method)
		case 'char'
			if ~any(strcmp(METHOD_ALL, method))
				error('%s not found in possible methods: %s', method, METHOD_ALL{:})
			end
		case 'cell'
			for m = method(:)
				if ~any(strcmp(METHOD_ALL, m))
					error('%s not found in possible methods: %s', method, METHOD_ALL{:})
				end
			end
		otherwise
			error('Second input (method) should be of type ''char'' or ''cell''. Use [] to use default method.')
	end
		
end

if ~exist('PLOT', 'var')
	PLOT = true;
end

%% Compute WPLI
if any(strcmpi(method, 'wpli'))
	
	%% Downsample lfp
	DOWNSAMPLED_FREQ = 1e3;  % Hz
	ds_step = floor(mea.SamplingRate / mea.skipfactor / DOWNSAMPLED_FREQ);
	DOWNSAMPLED_FREQ = mea.SamplingRate / mea.skipfactor / ds_step;
	lfp = mea.lfp;
	lfp = lfp(1:ds_step:end, :);
	Time = mea.Time();
	Time = Time(1:mea.skipfactor*ds_step:end);

	%% Set parameters for coherence calculations
	W = 5;  % bandwidth (Hz)
	T = 1;  % time of sample (s)
	p = 1;  % ntapers = W*T - p


	params.Fs = DOWNSAMPLED_FREQ;
	params.tapers = [W T p];
	params.pad = 0;  % padding factor for the FFT. if PAD = 0, we pad the FFT to the next power of 2
	params.fpass = [2 30];  % frequency band to be used in the calculation 
	% params.err = [1 .005];


	%% Compute the multitaper complex cross spectrum at 160 ms intervals
	window = 1;  % length of window (s)
	winSamples = window * DOWNSAMPLED_FREQ;  % samples per window
	winstep = floor(160e-3 * DOWNSAMPLED_FREQ);  % step through the data in 160 ms increments
	nSamples = size(lfp, 1) - winSamples + 1;  % number of samples for which to compute coherence

	% initialize vector to store coherence
	wpli_mean = zeros(size(1:winstep:nSamples));  

	i = 0;  % increment by one
	for ti = 1:winstep:nSamples  % step through time at winstep increments
		i = i + 1;

		% Compute the complex cross spectrum and frequencies
		[Sc,~,~,~,~,freq] = CrossSpecMatc(...
			lfp(ti : ti + winSamples-1, :), window, params);

		% Compute the weighted phase locking index for each channel pair
		wpli = squeeze(abs(mean(imag(Sc), 1)) ./ mean(abs(imag(Sc)), 1));

		% take the mean of the squared wpli values
		wpli_mean(i) = mean(wpli(:).^2, 'omitnan');
	end
	
	% Save results
	disp('Saving WPLI')
	mea.wpli = wpli_mean;
	mea.wpliFreq = freq;
	output.wpli = wpli_mean;
	output.wpliFreq = freq;

	if PLOT
		figure(); clf; fullwidth()
		plot(Time(1:winstep:ti), wpli_mean);
		title([strrep(mea.Name, '_', ' '), ' Coherence'])
		ylabel('Mean square WPLI')
		xlabel('Time (s)')
		axis tight
	end
end
	

%% Compute Fano factor
if any(strcmpi(method, 'fano'))
	output = mua_fano_factor(mea);
end

