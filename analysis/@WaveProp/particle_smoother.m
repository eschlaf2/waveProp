function [ax, H, T] = ...
        particle_smoother(W, ax, NPARTS, SPREAD_FACTOR, SD_WIN, ANG_RES, ...
        SMOOTHING_WIN, SHOW_DIRECTIONS, SHOW_EVOLUTION)
    


if nargin < 2 || isempty(ax), ax = axes(figure); end
if nargin < 3 || isempty(NPARTS), NPARTS = 1e4; end
if nargin < 4 || isempty(SPREAD_FACTOR), SPREAD_FACTOR = 1; end
if nargin < 5 || isempty(SD_WIN), SD_WIN = 100;  end % [in number of 10ms increments; i.e. 100 means use 1s of surrounding time]
if nargin < 6 || isempty(ANG_RES), ANG_RES = 128; end
if nargin < 7 || isempty(SMOOTHING_WIN), SMOOTHING_WIN = 1; end
if nargin < 8 || isempty(SHOW_DIRECTIONS), SHOW_DIRECTIONS = false; end
if nargin < 9 || isempty(SHOW_EVOLUTION), SHOW_EVOLUTION = false; end


Vx = W.Vx;
Vy = W.Vy;
V = [Vx(:) Vy(:)];
T = W.time;
N = numel(T) - 1;
% inds = find(T >= 60 & T <= 65);
inds = find(T >= -Inf & T <= Inf);  inds = inds(1:end-1);

NSTEPS = min(numel(inds), N);
sts = ax.NextPlot;


% mag = vecnorm([Vx' Vy']);
% sd = nanstd(mag);


dV = diff(V);
sd = fillmissing(movstd(dV, SD_WIN, 'omitnan'), 'nearest');
sd(sd == 0) = nan;
sd = fillmissing(sd, 'constant', nanstd(dV));

dDist = vecnorm(dV, 2, 2);   % 2-norm along 2nd dimension
sd_dist = fillmissing(movstd(dDist, SD_WIN, 'omitnan'), 'nearest');
sd_dist(sd_dist == 0) = nan;
sd_dist = fillmissing(sd_dist, 'constant', nanstd(dDist));


dphi = diff_phase(atan2(Vy, Vx));

% p_phi =@(ii, diffs) ksdensity(dphi(inds(ii) - SD_WIN:inds(ii)), diffs, 'bandwidth', pi/128); 
[p_phi, x] = ksdensity(dphi(~isnan(dphi)), 'bandwidth', pi/128, 'npoints', 1e3); 
dphi = fillmissing(dphi, 'nearest');
% sd_phi =@(ii) circ_std(dphi(inds(ii)-SD_WIN:inds(ii)));



% dV(isnan(dV(:, 1)), :) = [];
% dZ = diff_phase(Z);

first_ind = find(~isnan(Vx(inds)), 1);
samples = [Vx(inds(first_ind)), Vy(inds(first_ind))] ...
    + randn(NPARTS, 2) .* sd(1, :) * SPREAD_FACTOR;
S = nan(NSTEPS, NPARTS, 2);
v_est = nan(NSTEPS, 2);
phi_est = nan(NSTEPS, 3);
P = nan(NSTEPS, 1);


if SHOW_EVOLUTION
    ax.NextPlot = 'replacechildren';
    xlim(quantile(Vx(inds), [.05 .95]));
    ylim(quantile(Vy(inds), [.05 .95]));
    sc_samples = scatter(ax, samples(:, 1), samples(:, 2), 'filled');
    ax.NextPlot = 'add';
    sc_obs = scatter(ax, Vx(1), Vy(1), 'filled');
    sc_v_est = scatter(ax, Vx(1), Vy(1), [], [0 0 0], ...
        'marker', 'x', 'linewidth', 4);
    ax.NextPlot = 'replacechildren';
end

wrn = warning;
warning('off');
REDISTRIBUTE = 0;
for ii = first_ind:NSTEPS
    
    t_win = inds(ii) + (1:SD_WIN) - floor(SD_WIN/2);
    if t_win(1) < 1, t_win = t_win - t_win(1) + 1; end
    if t_win(end) > N, t_win = t_win - (t_win(end) - N); end
    
    
    PLOT = SHOW_EVOLUTION && logical(~mod(ii, 10));
    if ~mod(ii, 100)  % update p_phi every second
        [p_phi, x] = ksdensity(dphi(t_win), ...
            'npoints', 1e3, 'bandwidth', pi/128);
    end
    
    
    obs = [Vx(inds(ii)), Vy(inds(ii))];
    if REDISTRIBUTE  % redistribute particles
        dat = V(t_win, :);
        dat(isnan(dat(:, 1)), :) = [];
        samples = dat(randi(size(dat, 1), NPARTS, 1), :);
        REDISTRIBUTE = 0;
        fprintf('%d, ', ii)
    end
    
    samples = samples + randn(size(samples)) .* sd(ii, :)*SPREAD_FACTOR;

    if isnan(obs(1))
        S(ii, :, :) = samples;
        continue
    end
    
    dist = obs - samples;
    dist = sqrt(sum(dist.^2, 2));
    p_dist = normpdf(dist, 0, sd_dist(ii));


    diffs = diff_phase([repmat(atan2(obs(2), obs(1)), NPARTS, 1), ...
        atan2(samples(:, 2), samples(:, 1))], [], 2);
%         p = normpdf(diffs, 0, sd_phi(ii)) + 1e-6;
%         p = p_phi(ii, diffs) + 1e-6;
    p_angle = interp1(x, p_phi, diffs, 'nearest', 'extrap');
    p = p_dist + p_angle + 1e-6;
    P(ii) = nanmax(p) / 2;
    if P(ii) < .1; REDISTRIBUTE = true; end
    p = p./sum(p);  


    % resample
    new_inds = interp1(cumsum(p), (1:NPARTS)', rand(NPARTS, 1), ...
        'nearest', 'extrap');
    samples = samples(new_inds, :);

    if PLOT
        sc_samples.XData = samples(:, 1);
        sc_samples.YData = samples(:, 2);
        sc_samples.CData = p(new_inds);
        sc_obs.XData = obs(1);
        sc_obs.YData = obs(2);
        title(ax, num2str(ii))
        drawnow
    end


    
    
    S(ii, :, :) = samples;

    
    % compute estimates
%     v_est(ii, :) = sum(samples .* p(new_inds)./sum(p(new_inds)));
    [~, max_ind] = sort(p(new_inds), 'descend');
    v_est(ii, :) = mean(samples(max_ind(1:NPARTS/100), :));
    
    phi_samp = atan2(samples(:, 2), samples(:, 1));
    [mn, ul, ll] = circ_mean(phi_samp, p(new_inds));
    phi_est(ii, :) = [mn, ul, ll];
%     phi_ci(ii, 1) = 
    
    if PLOT
        sc_v_est.XData = v_est(ii, 1);
        sc_v_est.YData = v_est(ii, 2);
    end
end
warning(wrn);
ax.NextPlot = sts;


%% Visualize results
phi_temp = atan2(S(:, :, 2), S(:, :, 1));
rotate_by = W.RotateBy(1);
if rotate_by, phi_temp = angle(exp(1j*(phi_temp + rotate_by))); end
% phi2 = atan2(v_est(:, 2), v_est(:, 1));
% phi2(isnan(phi_est(:, 2))) = nan;
edges = linspace(-pi, pi, ANG_RES + 1);
H = arrayfun(@(ii) ...
    histcounts(phi_temp(ii, :), edges), ...
    1:NSTEPS, 'uni', 0);
H = cat(1, H{:}) / NPARTS;
Hlo = resample(H, 1, 10);
if SMOOTHING_WIN
    Hlo = smoothdata(Hlo, 1, 'gaussian', SMOOTHING_WIN, ...
        'samplepoints', T(inds(1:10:end))); 
end
aa = edges(1:end-1) + pi/ANG_RES;
contourf(ax, T(inds(1:10:end)), aa, Hlo' / (aa(2) - aa(1)), ...
    'levellist', linspace(.25, .75, 10));
if SHOW_DIRECTIONS
    
    ax.NextPlot = 'add';
%     plot(T(inds), phi_est(:, 1), '.');
    plot(T, W.Z, '.');
    ax.NextPlot = sts;
end
colormap(ax, 1-gray);

T = T(inds);


end