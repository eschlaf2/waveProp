function mea = exclude_channels(mea)
	% Exclude channels with mean firing rate less than 6 spikes per second 
	% during seizure or that do not show significant increases in firing rate 
	% during the first half of the seizure

	try firingRate = mea.firingRate; catch, firingRate = mua_firing_rate(mea); end
	time = mea.Time();
	seizure_inds = logical((time > 0) .* (time < time(end) - mea.Padding(2))); % Indicate when seizure is occurring
	exclude_channels = find(mean(firingRate(seizure_inds, :)) < 6);            % Only use channels with mean firing rate at least 6 spikes per second
	skipfactor = round(mea.SamplingRate / 100);                                % downsample to ~100 Hz
	fr = smoothdata(downsample(firingRate, skipfactor), 1, 'movmean', 1e3);    % smooth over 10 s windows
	time = downsample(time, skipfactor);  % downsample time to match fr

	mn = mean(fr(time < 0, :));    % get mean firing rate at baseline
	sd = std(fr);                  % ... and standard deviation over full recording
	fr = (fr - mn) ./ sd;          % normalize smoothed firing rate

	% Find channels that don't have significant increases in firing rate during 
	% the first half of the seizure
	halfWay = (time - mea.Padding(2)) / 2;
	exclude_channels = unique([exclude_channels find(max(fr(time < halfWay, :)) < 2, 1)])';

	% Add the channels to the list of bad channels
	ch = 1:size(mea.Data, 2);                                            % All channels
	ch(mea.BadChannels) = [];                                            % ... minus those already excluded
	mea.BadChannels = sort([mea.BadChannels; ch(exclude_channels)']);  % ... gives all bad channels
	disp('Channels excluded from analysis:')
	disp(ch(exclude_channels));

end
