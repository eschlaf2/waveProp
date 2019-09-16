k = 1000;
t0 = find(mea.Time >= 2, 1);
params.pad = 1;
params.tapers = [3 5];
[samples, trials, ~] = size(data);
direction = 'rows'; 

r = nextpow2(samples-1) + params.pad;
[C, phi, S12, S1, s2, phistd] = deal(nan(2^r+1,8,k));

for ii = 1:k
	if rem(ii, 1e4) == 0
		fprintf('ii=%d/%d\n', ii, k)
	end
	sub = squeeze(data(:, :, t0+ii));  % - mean(data(:, :, t0 + (1:k)), 3);;
	switch direction
		case 'columns'
			% columns
			ref = sub(:, 1) * ones(1, 8); 
			test = sub(:, 2:end);
		case 'rows'
			% rows
			ref = sub(1, :)' * ones(1, 8);
			test = sub(2:end, :)';
	end
	
	[C(:, :, ii), phi(:, :, ii), S12(:, :, ii), S1(:, :, ii), ...
		s2(:, :, ii), f, confC, phistd(:, :, ii)] = ...
		coherencyc(ref, test, params);
end
%%	
pcolor = @emilys_pcolor;
unwrapZ = @(X) unwrap(2.2*X)/2.2;
% unwrapZ = @(X) unwrap(1*X)/1;
diffS = @(X) [mean(diff(X)); diff(X)];

fignum = 5 + strcmpi(direction, 'rows');
figure(fignum); fullwidth(); 
subplot(121); 
pcolor(1:8, f, mean(C, 3, 'omitnan')', 'clims', [confC(1), 1]);
subplot(122); 
phi(C < confC(1)) = nan;
phi = unwrap(phi, [], 3);
res = unwrapZ(mean(phi, 3, 'omitnan'))./(1:8);
% res = unwrapZ(sum((phi .* C)./sum(C, 3), 3, 'omitnan'));
% res = median(phi, 3, 'omitnan');

% pcolor(1:8, f, mod(res, 2*pi)', 'cmap', 1-hsv, 'clims', [0 2*pi]);
pcolor(1:8, f, res')
% pcolor(1:8, f, diffS(res)'/mean(diff(f)));
l = legend('location', 'southwest'); l.Visible = 'off';
axes('Position', l.Position .* [1.02 1.2 1 2])
plot(f, mean(res, 2, 'omitnan'));
% ylim([-pi pi])

%%

figure(4); 
for ii = 1:10:k 
	sub = data(:, :, t0+ii);  % - mean(data(:, :, t0 + (1:k)), 3);
	pcolor(sub); 
	title(num2str(t0+ii)); 
	drawnow(); pause(.1); 
end

%%

% phi(C < confC(1)) = nan; 
% phi = unwrap(phi, [], 3);
% delays.rows.phi = phi;
% delays.rows.C = C;
% delays.cols.phi = 


X = squeeze(mean((smoothdata(delays.cols.phi, 3, 'movmean', 200)./(1:8) ), 2, 'omitnan'));
Y = squeeze(mean((smoothdata(delays.rows.phi, 3, 'movmean', 200)./(1:8) ), 2, 'omitnan'));
dir = atan2(Y, X);
mag = sqrt(X.^2 + Y.^2);
mask = mag > 1 | mag == 0;
mag(mask) = nan; 
dir(mask) = nan;

[m, mx] = sort(std(diff(unwrap(dir'))));
finds = mx(1:5);

figure(9); clf
pcolor(mea.Time, f(finds), dir(finds, :)', 'cmap', hsv, 'clims', [-pi pi])
% yyaxis('right'); 
% colormap(hsv)
% trend = smoothdata(median(dir(finds, :), 'omitnan'), 'movmed', 100);
% set(gca, 'colormap', hsv, 'nextplot', 'replacechildren')
% h = scatter(mea.Time, trend, 20, trend);
% h.Parent.CLim = [-pi pi];

figure(10); clf
pcolor(mea.Time, f, log(mag'), 'clims', [-10 1])
% yyaxis('right'); 
% plot(mea.Time, mean(mag, 'omitnan'))
% % pcolor(mea.Time, f, atan(Y./X)', 'cmap', hsv)

