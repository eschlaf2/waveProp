function mea = mua_fano_factor(mea, PLOT, method)
% From Schevon et al., 2012 - details not published, but that the FF
% decreases to below baseline while the firing rate is still high once the
% region is recruited

if ~exist('PLOT', 'var')
	PLOT = true;
end

if ~exist('method', 'var')
	method = 2;
end

switch method
	case 1  % ff of sum of all events (fine grain - step = 1)
%%
		windMS = 100;  % ms window to use
		window = mea.SamplingRate * 1e-3 * windMS;  % samples per ms * ms to use

		fr = padarray(sum(mea.events, 2), window, 'symmetric');

		ff = arrayfun(@(i) var(fr(i - window : i + window - 1)) / ...
			mean(fr(i - window + 1 : i + window)), ...
			(1:size(mea.firingRate, 1)) + window);
%%	
	case 2  %% mean ff of firing rate on each channel (course grain: step = window)
%%
		windMS = 1e3;  % ms window to use
		window = mea.SamplingRate * 1e-3 * windMS;  % samples per ms * ms to use
		[numSamples, numCh] = size(mea.firingRate);  % array size

% 		fr = sum(mea.events, 2);
% 		fr2 = reshape(fr(1:floor(numel(fr) / window) * window), window, []);
		fr2 = reshape((mea.firingRate(1:(floor(numSamples / window) * window), :)), ...
			window, [], numCh);
		ff2 = squeeze((var(fr2) + 1) ./ mean(fr2 + 1));
		
		if PLOT
			figure(3); plot(mea.Time(window:window:end), mean(ff2, 2))
			title([strrep(mea.Name, '_', ' ') ': Fano factor'])
			xlabel('Time (s)');
			ylabel('FF')
			axis('tight')
		end
%%
	case 3  % ff computed for events on each channel and then averaged (course grain)
%%
		windMS = 1e3;  % ms window to use
		window = mea.SamplingRate * 1e-3 * windMS;  % (samples per ms) * (ms to use)
		[numSamples, numCh] = size(mea.firingRate);  % array size

		fr = mea.events(1:(floor(numSamples / window) * window), :);  % prep to reshape
		fr = reshape(fr, window, [], numCh);  % reshape array into time bins

		vr = squeeze(var(fr));  % get sdev of each bin
		mn = squeeze(mean(fr + 1));  % ... and mean of each bin

		ff = vr ./ mn;  % compute fano factor at each time point for each channel
		
		if PLOT
			figure(3);
			plot(mea.Time(window:window:end), mean(ff, 2, 'omitnan'))  % plot the mean ff
		end
%%		
end

