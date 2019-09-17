function [] = group_velocity(mea, wavetimes, halfwin)
% Input: wavetimes and halfwin in ms

% k = 1000;
% t0 = find(mea.Time >= 2, 1);
nW = numel(wavetimes);
time = mea.Time() * 1e3;  % time in ms
k = 2 * round(halfwin / mea.SamplingRate * 1e3) + 1;
position = mea.Position;
data = data2grid(mea.Data, position);
ddims = diff(size(data));
if ddims(1) < 0, data = [data data(:, 1:-ddims(1), :) * 0]; end
if ddims(1) > 0, data = [data data(1:ddims(1), :, :) * 0]; end
params.pad = 1;
params.tapers = [3 5];
samples = size(data, 1);
% samples = size(data, 1);
% direction = {'rows', 'cols'}; 

	r = nextpow2(samples-1) + params.pad;
	[C, phi, S12, S1, s2, phistd] = deal(nan(2^r+1,8,k,nW));
for w = 1:length(wavetimes)
	t0 = find(time >= wavetimes(w), 1);
	for ii = -halfwin:halfwin
		if rem(ii, 1e4) == 0
			fprintf('ii=%d/%d\n', ii, k)
		end
		sub = squeeze(data(:, :, t0+ii));  % - mean(data(:, :, t0 + (1:k)), 3);;
% 		switch d
% 			case 'columns'
				% columns
				refc = sub(:, 1) * ones(1, trials-1); 
				testc = sub(:, 2:end);
% 			case 'rows'
				% rows
				refr = sub(1, :)' * ones(1, trials-1);
				testr = sub(2:end, :)';
% 		end

		[C(:, :, ii, w), phi(:, :, ii, w), S12(:, :, ii, w), S1(:, :, ii, w), ...
			s2(:, :, ii, w), f, confC, phistd(:, :, ii, w)] = ...
			coherencyc([refc refr], [testc testr], params);
	end
end

%%
phi = unwrap(phi, [], 3);  % unwrap phi 
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

