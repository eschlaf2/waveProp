% function preproc_data(mea)

data = double(mea.Data);
[nT, nCh] = size(data);
[coeff, score, ~, ~, explained, ~] = pca(data);

%% Figure 10: raw v. preprocessed
figure(10); clf; fullwidth(1)
[r, c] = deal(2);

subplot(r, c, 1)  % Variance explained by PCs in raw data
stem(cumsum(explained(1:20)));
hold on;
stem(explained(1:20))
hold off;
title('Variance explained (Raw)')
xlabel('Component'); ylabel('Percent explained')

subplot(r, c, 2)  % Standard deviation of raw traces
stem(std(data))
title('Standard deviation (Raw)')
xlabel('Channel'); ylabel('STD');

%% Preprocessing
% Use PCA to separate neural signal from electrode drift(?) 
hiVar = find(abs(zscore(std(data))) > 1);  % Channels with high variance drown out neural signal
[mCh, singleCh] = find(abs(coeff) > .6);  % Find PCs that describe the activity of a single channel
badPCs = singleCh(any(hiVar' == mCh'));  % Filter to channels with high variance and big PC contributions (large changes in electrode reading that are probably not biological)
badPCs = union(badPCs, find(explained > 60));
PCs = 1:nCh; PCs(badPCs) = [];  % Exclude these PCs
data = score(:, PCs) * coeff(:, PCs)';  % ... and rebuild the data

%% Figure 13: Scores
figure(13); clf; fullwidth(true)
cmap = lines;
goodPCs = 1:nCh; goodPCs(badPCs) = [];
plot(zscore(score(:, goodPCs))/10 + goodPCs, 'color', cmap(1, :)); hold on;
plot(zscore(score(:, badPCs))/10 + badPCs', 'color', cmap(2, :)); hold off;
ylim([0 nCh+1])

%% Rerun PCA
[coeff, score, latent, ~, explained, ~] = pca(data);  % Rerun PCA

%% Figure 11: PCs of preprocessed data
figure(11); clf; fullwidth(1);
for ii = 1:19
	subplot(4, 5, ii)
	imagesc(sparse(mea.Position(:, 1), mea.Position(:, 2), coeff(:, ii)))
	axis square
	title(sprintf('PC%d: %0.2f%%', ii, explained(ii)));
end
subplot(4, 5, 20);
colorbar('Location', 'southoutside')
set(gca, 'Visible', 'off')
annotation(figure(11),'textbox',...
    [0.843 0.133 0 0.096],...
    'String',strrep(mea.Name, '_', ' '),...
	'LineStyle', 'none', ...
    'FontWeight','bold',...
    'FontSize',18,...
	'HorizontalAlignment', 'center', ...
    'FitBoxToText','on');

%% Figure 10: Raw v. pre-processed (cont)
figure(10);
subplot(r, c, 3)  % Variance explained by each PC (preproc'd data)
stem(cumsum(explained(1:20)));
hold on;
stem(explained(1:20))
hold off;
title('Variance explained (Preproc)')
xlabel('Component'); ylabel('Percent explained')

subplot(r, c, 4)  % STD on each channel (preproc'd data)
stem(std(data))
title('Standard deviation (Preproc)')
xlabel('Channel'); ylabel('STD');

%% Figure 12: Remove channels with extreme deviations
figure(12); clf; fullwidth(true)
cmap = lines;
badCh = union(find(abs(zscore(std(data))) > 2), ...  % Channels with extreme variance should be excluded
	find(std(data) < 1));  % ... as well as inactive channels
goodCh = 1:nCh; goodCh(badCh) = [];
plot(zscore(data(:, goodCh))/10 + goodCh, 'Color', cmap(1, :)); hold on;  % Plot good channels
plot(zscore(data(:, badCh))/10 + badCh, 'Color', cmap(2, :)); hold off  % ... and bad channels
title({'Bad Channels:'; num2str(badCh)})
axis('tight')
ylim([0 nCh+1])

%% Save results
mea.BadChannels = badCh(:);
save(mea.Path, '-v7.3', '-struct', 'mea')

print(10, [mea.Name '_raw_v_preproc'], '-dpng')
print(11, [mea.Name '_PCs_final'], '-dpng')
print(12, [mea.Name '_preproc'], '-dpng')
print(13, [mea.Name '_scores'], '-dpng')



