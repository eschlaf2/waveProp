function [] = plot_Neuroport(Data, samplingRate)

[nSamp, nCh] = size(Data);
	if nCh > 96
		Data = Data(:, 1:96);
		nCh = 96;
	end
% 	Data = Data ./ max(abs(Data)) / 2 + int16(1 : nCh);
% 	Data = smoothdata(Data, 'gaussian', samplingRate / 10);
	spikes = abs(zscore(Data)) > 4;
% 	Data = Data ./ max(abs(Data));
	T = (1:nSamp) / samplingRate;
	figure(998); fullwidth(); 
% 	plot((1:nSamp-1) / samplingRate, zscore(diff((Data))) + (1:nCh))
	yyaxis right
	plot(T, smoothdata(std(Data, [], 2), 'movmean', samplingRate), 'k', 'linewidth', 2); 
	yyaxis left
	set(gca, 'colororder', cool(nCh)); set(gca, 'nextplot', 'replacechildren')
	plot(T, (spikes) + (1:nCh))
	axis tight;