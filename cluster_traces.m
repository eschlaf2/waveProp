% Cluster traces in data (TxCh) variable using kmeans (with k=15).

kN = 15;
[idx, C] = kmeans(zscore(data)', kN);
figure(); 
subplot(2,2,1); plot(C'/10 + (1:kN));

dist = pdist(C);
subplot(2,2,2); imagesc(squareform(dist));
subplot(223); stem(sum(squareform(dist)))
subplot(224); colormap(hsv(kN + 1));
imagesc(sparse(position(:, 1), position(:, 2), idx));