function [mea] = mua_events(mea)

intervalM = mea.SamplingRate / 1e3;  % samples per ms
mn = mean(mea.mua(mea.Time < 0, :), 1);
sd = std(mea.mua(mea.Time < 0, :), [], 1);
temp = (mea.mua - mn) ./ sd;
mea.artefacts = find(abs(temp) > 16);
mea.mua(mea.artefacts) = 0;
temp(temp > 0) = 0;
mea.events = false(size(temp));
for ch = 1:size(temp, 2)
	temp = mea.mua(:, ch);
	[~, inds] = findpeaks(-zscore(temp), ...
		'minpeakdistance', intervalM, 'minpeakheight', 4);  
	mea.events(inds, ch) = true;
end