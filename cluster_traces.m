% Cluster traces in data (TxCh) variable using kmeans (with k=15).

kN = 15;
[idx, C] = kmeans(zscore(data)', kN);
figure(); plot(C'/10 + (1:kN));

dist = pdist(C);
figure(); imagesc(squareform(dist));
figure(); stem(sum(squareform(dist)))