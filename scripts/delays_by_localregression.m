% delays_by_localregression
% basename = compute_coherograms(pat, seizure, T, W, DS, units, toi);  %
% (See quickscript)
load(basename, ...  % Load the following from basename
    'f', 'params', 'phi', 'C', 't', 'confC', 'pat', ...
    'seizure', 'position', 'pairs', 'units')

df = mean(diff(f)); 
winsz = 3 * units;  % Hz
thresh = 5e-2; 
fband = params.fpass;
if ~exist('toi', 'var'); toi = [-Inf Inf]; end
if ~exist('tau', 'var'); tau = 'group'; end  % 'group' or 'phase'
% MASK = false;

if isinteger(phi); phi = -single(phi) / 1e4; end
finds = (f >= fband(1)) & (f <= fband(2));
tinds = (t >= toi(1)) & (t <= toi(2));
t = t(tinds); f = f(finds);

[nt, nf, np] = size(C(tinds, finds, :));

transform = @(A) reshape(permute(A(tinds, finds, :), [2 1 3]), nf, []);
Cf = transform(C);  % limit to band of interest
phif = transform(phi);  % ... same for phi
clear C phi

%% Delays
disp('Computing delays')

delays = dblr(phif, f, tau, Cf, confC, winsz);

if exist('MASK', 'var') && MASK  % Keep only the longest streak of significant data points
	delays = keep_streak(delays, winsz / df);
end

delaysR = reshape(delays, nf, nt, []);
clear delays phif


%% All wave directions
disp('Computing directions')

MIN_RATIO_FINITE = 0.2;
pos = position(pairs(:, 2), :);

nslots = str2double(getenv('NSLOTS'));  % Check for parpool
if ~isnan(nslots), p = parpool(nslots); end

warning('off', 'stats:statrobustfit:IterationLimit');
H = [0 1 0; 0 0 1];  
c = [0; 0];

if isnan(nslots)
    disp('Using arrayfun ...')
    tic %#ok<*UNRCH>
    delaysR2 = reshape(delaysR, [], np)';  % reshape the delays so that each column is a single time-freq point and rows are pairs
    num_unique = arrayfun(@(ii) ...  % Number of unique, non-nan values in each column
        numel(unique(delaysR2(isfinite(delaysR2(:, ii)), ii))), 1:length(delaysR2));
    finite = sum(isfinite(delaysR2));  % number of non-nan, non-inf delays in each column
    Z = nan * num_unique;  % Initialize Z with nans
    inds = ...  % Need at least 3 unique values and MIN_RATIO_FINITE finite values
        find((finite > max(MIN_RATIO_FINITE * np, 3)) & (num_unique >= 3));  
    
    disp('Computing fits')
    tic
    [beta, stats] = ...  % Compute fits
        arrayfun(@(ii) robustfit(pos, delaysR2(:, ii), 'fair'), ...
        inds, 'uni', 0);  
    toc
    
    disp('Computing significance')
    tic
    pdel = ...  % ... and significance
        cellfun(@(A, B) linhyptest(A, B.covb, c, H, B.dfe), beta, stats);  
    toc
    
    indsB = find(pdel < thresh);  % Keep significant fits
    disp('Computing Z')
    tic
    Z(inds(indsB)) = ...  % Fill in wave directions for good fits
        arrayfun(@(ii) angle(pinv(beta{ii}([2 3])) * [1; 1i]), indsB);
    Z = reshape(Z, nf, nt);
    toc
    
    toc
else
    disp('Using parfor ...')
    tic
    [Z, pdel, pct, V] = deal(nan(nf, nt));
    parfor ii = 1:nf  % For each frequency
        if ~mod(ii, 100), fprintf('ii=%d/%d\n', ii, nf), end
        for jj = 1:nt  % ... and time point
            delays2fit = squeeze(delaysR(ii, jj, :));  % ... collect the delays for each pair
            num_unique = numel(unique(delays2fit(isfinite(delays2fit))));  % ... calculate the number of finite unique values
            if num_unique < 3, continue; end  % check that there are enough unique points to fit a plane
            finite = sum(isfinite(delays2fit));  % Calculate the number of finite delays
            pct(ii, jj) = finite / np;  % (FYI only)
            if finite <= max(MIN_RATIO_FINITE * np, 3); continue; end  % check enough delay data is not NaN.
            [beta,stats] = robustfit(pos, delays2fit, 'fair');  % fit the delay vs two-dimensional positions
            ptemp = linhyptest(beta, stats.covb, c, H, stats.dfe);  % Compute significance of fit
            if ptemp >= thresh || isnan(ptemp); continue; end  % Stop if fit is not significant
            Vt = pinv(beta(2:3));  % velocity is the pseudoinvervse of the fitted wave
            Z(ii, jj) = angle([1 1i] * Vt(:));  % Z is the angle of the velocity
            V(ii, jj) = abs([1 1i] * Vt(:));
            pdel(ii, jj) = ptemp;  % Store p-value
        end
    end
    toc
    [msg,msgID] = lastwarn
    delete(p)
end

%% Imagesc Z (angles computed using delays)
figure();
emilys_pcolor(t, f, Z', 'cmap', hsv(80), 'clim', [-pi,pi]);
line(t, 13 * ones(size(t)), 'color', 'black', 'linewidth', 2)
xlabel('Time (s)');
ylabel('Freq (Hz)')
title(sprintf('%s Seizure %d\n%s delays', pat, seizure, tau));

