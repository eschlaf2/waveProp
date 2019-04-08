nS = numel(wave_fit.phis);
X = mea.X; Y = mea.Y;
nCh = numel(X);
T = size(wave_fit.phis{1}, 3);

for t = 1:nS
	phis = wave_fit.phis{t};


	phis2 = nan(T, nCh);
	E = zeros(T, nCh);
	S = E;
	for i = 1:nCh
		phis2(:, i) = phis(X(i), Y(i), :);
		temp = diff([0; find(isnan(phis2(:, i))); T]);
		[E(i), ind] = max(temp - 1);
		S(i) = ind + sum(temp(1:ind-1) - 1);
		nan_mask = true(nCh, 1);
		nan_mask(S(i) + (0:E(i))) = false;
		phis2(nan_mask, i) = nan;
	end
	
	figure(1); plot(phis2/2/pi + (1:nCh))

% 	figure(5); clf; 
% 	imagesc(isnan(phis2)); axis xy; hold on; plot(E+S, 'r*'); plot(S, 'g*');
% 	title(num2str(t));
% 	pause();
end