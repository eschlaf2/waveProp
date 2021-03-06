function mea = exclude_channels(mea)
	% Exclude channels with mean firing rate less than 6 spikes per second 
	% during seizure or that do not show significant increases in firing rate 
	% during the first half of the seizure
	
	EXCLUDE = mea.params.exclude_rate;  % (default: 6 spikes/sec)
	REQUIRE_INCREASE = mea.params.require_increase;  % (default: true)
	
	try firingRate = mea.firingRate; catch, firingRate = mua_firing_rate(mea); end
	
	time = mea.Time();
	seizure_inds = logical((time > 0) .* (time < time(end) - mea.Padding(2))); % Indicate when seizure is occurring
	exclude_channels = find(mean(firingRate(seizure_inds, :)) < EXCLUDE);            % Only use channels with mean firing rate at least 6 spikes per second
	
	if REQUIRE_INCREASE
		% Find channels that don't have significant increases in firing rate during 
		skipfactor = round(mea.SamplingRate / 100);                                % downsample to ~100 Hz
		fr = smoothdata(downsample(firingRate, skipfactor), 1, 'movmean', 1e3);    % smooth over 10 s windows
		time = downsample(time, skipfactor);  % downsample time to match fr

		mn = mean(fr(time < 0, :));    % get mean firing rate at baseline
		sd = std(fr);                  % ... and standard deviation over full recording
		fr = (fr - mn) ./ sd;          % normalize smoothed firing rate

		% the first half of the seizure
		halfWay = (time - mea.Padding(2)) / 2;
		exclude_channels = unique([exclude_channels find(max(fr(time < halfWay, :)) < 2, 1)])';

	end
	% Add the channels to the list of bad channels
	ch = 1:size(mea.Data, 2);                                            % All channels
	ch(mea.BadChannels) = [];                                            % ... minus those already excluded
	mea.BadChannels = sort([mea.BadChannels(:); ch(exclude_channels)']);  % ... gives all bad channels
	if isfield(mea, 'lfp'), mea.lfp(:, exclude_channels) = []; end
	if isfield(mea, 'mua'), mea.mua(:, exclude_channels) = []; end
	disp('Channels excluded from analysis:')
	disp(ch(exclude_channels));

end
