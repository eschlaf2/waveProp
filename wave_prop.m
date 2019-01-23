%%
% All time units should be converted to ms for consistency


sig = .25;
halfWin = 50;  % half window around discharge event (ms)
plotTitle = 'BW09 Seizure 1';

%%

[~, waveTimes] = findpeaks(mean(mea.firingRate, 2), mea.SamplingRate / 1e3, ...
	'minpeakheight', 150, ...
	'minpeakprom', 50, 'minpeakdistance', 100);  % The frequency here is given in samples per ms
waveTimes = waveTimes - mea.Padding(1) * 1e3;  % Account for padding (in ms)
mea.TimeMs = mea.Time * 1000;
T = mea.TimeMs' .* mea.events;
T(T == 0) = nan;
numCh = numel(mea.X);
numWaves = numel(waveTimes);
% numWaves = 2;
beta = nan(3, numWaves);
V = nan(2, numWaves);
p = nan(1, numWaves);

%%
for i = 1:numWaves
	t = waveTimes(i);
	inds = find((mea.TimeMs >= t - halfWin) .* (mea.TimeMs <= t + halfWin));
	events = T(inds, :)';
	events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
	for ch = 1:numCh
		temp = events{ch};
		temp(isnan(temp)) = [];
		events{ch} = temp;
	end
	[beta(:, i), V(:, i), p(i)] = ...
		SpatialLinearRegression(events, [mea.X, mea.Y], ...
		'switch_plot', 0, 'Lossfun','L2');
end
Z = angle(complex(V(1, :), V(2, :)));
%%

figure(2); clf; fullwidth()
subplot(1, 4, 1:3);
plot(mea.TimeMs, zscore(mean(mea.firingRate, 2)), 'color', .5*[1 1 1]); 
hold on;
axis tight; ylim([-6, 6]); 
N = sqrt(sum(V.^2));
% Vn = beta(1:2, :);  
Vn = 2*V;
cmap = jet(numWaves); % cols(1, :) = [];
% plot([0; Vn(1, 1) + 0] + timeFr(waveTimes(1))/1e3, [0; Vn(2, 1)], ...
% 	'color', cmap((p(1) < sig) + 1, :)); 
grid on;
title(plotTitle)
xlabel('Time (ms)')

subplot(1, 4, 4); 
compass(V(1, p < sig), V(2, p < sig)); hold on;

for i = 2:numWaves
	if p(i) > sig
		continue
	end
	subplot(1, 4, 1:3);
% 	plot([0; Vn(1, i) * (size(firingRate, 1)/2e3)] + timeFr(waveTimes(i))/1e3, ...
% 		[0; Vn(2, i)], 'color', cmap(i, :), 'linewidth', 2); ylim([-6, 6]); % xlim([0 126]); 
	plot(waveTimes(i), Z(i), 'ko', 'MarkerFaceColor', cmap(i, :));
	
	subplot(1, 4, 4); 
	h = compass(V(1, i), V(2, i)); hold on
	h.Color = cmap(i-1, :);
	h.LineWidth = 2;
	
	drawnow()
	pause(1e-2)
end
hold off
