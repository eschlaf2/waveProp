%%

sig = .05;
plotTitle = 'BW09 Seizure 1';

%%

[~, waveTimes] = findpeaks(smooth(mean(firingRate,2)), 'minpeakheight', 150, 'minpeakprom', 50, 'minpeakdistance', 4);
mea.TimeMs = mea.Time * 1000;
timeFr = linspace(mea.TimeMs(1), mea.TimeMs(end), size(firingRate, 1));
T = mea.TimeMs' .* mea.events;
T(T == 0) = nan;
numCh = numel(mea.X);
numWaves = numel(waveTimes);
beta = nan(3, numWaves);
V = nan(2, numWaves);
p = nan(1, numWaves);

%%
for i = 1:numWaves
	t = waveTimes(i);
	inds = find((mea.TimeMs >= t - 100) .* (mea.TimeMs <= t + 100));
	events = T(inds, :)';
	events = mat2cell(events, ones(size(events, 1), 1), size(events, 2));
	for ch = 1:numCh
		temp = events{ch};
		temp(isnan(temp)) = [];
		events{ch} = temp;
	end
	[beta(:, i), V(:, i), p(i)] = SpatialLinearRegression(events, [mea.X, mea.Y], 'switch_plot', 0, 'Lossfun','L2');
end

%%

figure(2); clf; fullwidth()
subplot(1, 4, 1:3);
plot(timeFr / 1e3, zscore(mean(firingRate, 2)), 'color', .5*[1 1 1]); hold on;
axis tight; ylim([-6, 6]); 
N = sqrt(sum(V.^2));
% Vn = beta(1:2, :);  
Vn = 2*V;
cmap = jet(125); % cols(1, :) = [];
plot([0; Vn(1, 1) + 0] + timeFr(waveTimes(1))/1e3, [0; Vn(2, 1)], 'color', cmap((p(1) < sig) + 1, :)); hold on;
grid on;
title(plotTitle)
xlabel('Time (s)')

subplot(1, 4, 4); 
compass(V(1, p < sig), V(2, p < sig)); hold on;

for i = 2:numWaves
	if p(i) > sig
		continue
	end
	subplot(1, 4, 1:3);
	plot([0; Vn(1, i) * (size(firingRate, 1)/2e3)] + timeFr(waveTimes(i))/1e3, ...
		[0; Vn(2, i)], 'color', cmap(i, :), 'linewidth', 2); ylim([-6, 6]); % xlim([0 126]); 
	plot(timeFr(waveTimes(i))/1e3, Z(i), 'ko', 'MarkerFaceColor', cmap(i, :));
	
	subplot(1, 4, 4); 
	h = compass(V(1, i), V(2, i)); hold on
	h.Color = cmap(i-1, :);
	h.LineWidth = 2;
	
	drawnow()
	pause(1e-2)
end
hold off
