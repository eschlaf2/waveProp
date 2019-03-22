function [] = plot_Neuroport(Data, samplingRate)
% Plots large deviations in Data (more than 4 sd from mean). Traces are
% stacked and the standard deviation across the electrodes at each point in
% times is shown in a thick black line.

[nSamp, nCh] = size(Data);
if nCh > 96
	Data = Data(:, 1:96);
	nCh = 96;
end
spikes = abs(zscore(Data)) > 4;
T = (1:nSamp) / samplingRate;
figure(998); fullwidth(); 
yyaxis right
plot(T, smoothdata(std(Data, [], 2), 'movmean', samplingRate), 'k', 'linewidth', 2); 
yyaxis left
set(gca, 'colororder', cool(nCh)); set(gca, 'nextplot', 'replacechildren')
plot(T, (spikes) + (1:nCh))
axis tight;