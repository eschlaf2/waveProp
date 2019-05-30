% function preproc_data(mea)

raw = double(mea.Data);
[nT, nCh] = size(raw);
badCh = find(arrayfun(@(ii) numel(unique(raw(:, ii))), 1:nCh) < 164);  % Values are in range +/- 8192; 164 ~ 1% of range
goodCh = 1:nCh; goodCh(badCh) = [];
[coeff, score, ~, ~, explained, ~] = pca(raw(:, goodCh));
nPCs = size(coeff, 2);

%% Preprocessing
% Use PCA to separate neural signal from electrode drift(?) 
hiVar = find(abs(zscore(std(raw))) > 1);  % Channels with high variance drown out neural signal
[mCh, singleCh] = find(abs(coeff) > .6);  % Find PCs that describe the activity of a single channel
badPCs = singleCh(any(hiVar' == goodCh(mCh)));  % Filter to channels with high variance and big PC contributions (large changes in electrode reading that are probably not biological)
% badPCs = union(badPCs, find(explained > 60));
PCs = 1:nPCs; PCs(badPCs) = [];  % Exclude these PCs
data = raw;
data(:, goodCh) = score(:, PCs) * coeff(:, PCs)';  % ... and rebuild the data

%% Figure 13: Scores
figure(13); clf; fullwidth(true)
cmap = lines;
plot(zscore(score(:, PCs))/10 + PCs, 'color', cmap(1, :)); hold on;
plot(zscore(score(:, badPCs))/10 + badPCs', 'color', cmap(2, :)); hold off;
axis tight;
xlabel('Time (samples)'); ylabel('PC');
ylim([0 nCh+1])
title({[strrep(mea.Name, '_', ' ') ' Bad PCs:']; num2str(badPCs)})

%% Figure 11: Show first 20 PCs
figure(11); clf; fullwidth(1);
for ii = 1:19
	subplot(4, 5, ii)
	imagesc(sparse(mea.Position(goodCh, 1), mea.Position(goodCh, 2), coeff(:, ii)))
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

%% Rerun PCA
% [coeff, score, latent, ~, explained, ~] = pca(data(:, goodCh));  % Rerun PCA
explainedPreproc = explained(PCs);
explainedPreproc = explainedPreproc / sum(explainedPreproc);

%% Figure 10: Raw v. preprocessed

figure(10); clf; fullwidth(1)
[r, c] = deal(2);

subplot(r, c, 1)  % Variance explained by PCs in raw data
stem(cumsum(explained(1:20)));
hold on;
stem(explained(1:20))
hold off;
title([strrep(mea.Name, '_', ' ') 'Variance explained (Raw)'])
xlabel('Component'); ylabel('Percent explained')

subplot(r, c, 2)  % Standard deviation of raw traces
yyaxis left
stem(goodCh, std(raw(:, goodCh)))
ylabel('STD');
yyaxis right
stem(goodCh, zscore(std(raw(:, goodCh))));
ylabel('STD (normed)')
title('Standard deviation (Raw)')
xlabel('Channel'); 

figure(10);
subplot(r, c, 3)  % Variance explained by each PC (preproc'd data)
stem(cumsum(explainedPreproc(1:20)));
hold on;
stem(explainedPreproc(1:20))
hold off;
title('Variance explained (Preproc)')
xlabel('Component'); ylabel('Percent explained')

subplot(r, c, 4)  % STD on each channel (preproc'd data)
yyaxis left
stem(goodCh, std(data(:, goodCh)))
ylabel('STD');
yyaxis right
stem(goodCh, zscore(std(data(:, goodCh))));
ylabel('STD (normed)')
title('Standard deviation (Preproc)')
xlabel('Channel'); ylabel('STD');

%% Figure 12: Remove channels with extreme deviations
figure(12); clf; fullwidth(true)
cmap = lines;
badCh = union(goodCh(zscore(std(data(:, goodCh))) < -2), badCh);  % Channels with low variance should be excluded
goodCh = 1:nCh; goodCh(badCh) = [];
plot(zscore(double(mea.Data(:, goodCh))) ./ 10 + goodCh, 'Color', cmap(1, :)); hold on;  % Plot raw data
plot(zscore(data(:, goodCh))/10 + goodCh, 'Color', cmap(4, :)); hold on;  % ... and good channels
plot(zscore(data(:, badCh))/10 + badCh, 'Color', cmap(2, :)); hold off  % ... and bad channels
title({[strrep(mea.Name, '_', ' ') ' Bad Channels:']; num2str(badCh)})
axis('tight')
ylim([0 nCh+1])

%% Save results
% mea.BadChannels = badCh(:);
% save(mea.Path, '-v7.3', '-struct', 'mea')

fprintf('%d, ', badCh); fprintf('\b\b\n')

print(10, [mea.Name '_raw_v_preproc'], '-dpng')
print(11, [mea.Name '_PCs_final'], '-dpng')
print(12, [mea.Name '_preproc'], '-dpng')
print(13, [mea.Name '_scores'], '-dpng')

Name = mea.Name;
clear raw mea
save([Name '_preproc'])
