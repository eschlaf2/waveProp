function mea = mua_fano_factor(mea)
% From Schevon et al., 2012 - details not published, but that the FF
% decreases to below baseline while the firing rate is still high once the
% region is recruited

windMS = 100;  % ms window to use
window = mea.SamplingRate * 1e-3 * windMS;  % samples per ms * ms to use

fr = padarray(sum(mea.events, 2), window, 'symmetric');

% temp = zeros(2 * window + 1, size(mea.firingRate, 1));
% 
% for i = 1:(2 * window + 1)
% 	temp(i, :) = fr(i : end - (2 * window + 1) + i);
% end
% 
% ff = var(temp) ./ mean(temp);


ff = arrayfun(@(i) var(fr(i - window : i + window - 1)) / ...
	mean(fr(i - window + 1 : i + window)), ...
	(1:size(mea.firingRate, 1)) + window);

%%
windMS = 1e3;  % ms window to use
window = mea.SamplingRate * 1e-3 * windMS;  % samples per ms * ms to use

fr = sum(mea.events, 2);
fr2 = reshape(fr(1:floor(numel(fr) / window) * window), window, []);
ff2 = var(fr2) ./ mean(fr2);
figure(3); plot(mea.Time(1:window:end-window), ff2)

