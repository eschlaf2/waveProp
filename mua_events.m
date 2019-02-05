function [mea] = mua_events(mea)

intervalM = mea.SamplingRate / 1e3;  % samples per ms
mn = mean(mea.mua(mea.Time < 0, :), 1);  % get the mean of the baseline
sd = std(mea.mua);  % and the overall sd
temp = (mea.mua - mn) ./ sd;  % zscore based on baseline mean
mea.artefacts = (abs(temp) > 16);
mea.mua = mea.mua - (mea.mua - mn).* mea.artefacts;
% temp(temp > 0) = 0;
mea.events = false(size(temp));
for ch = 1:size(temp, 2)
	temp = mea.mua(:, ch);
	[~, inds] = findpeaks(-normalize(temp), ...
		'minpeakdistance', intervalM, 'minpeakheight', 4);  
	mea.events(inds, ch) = true;
end