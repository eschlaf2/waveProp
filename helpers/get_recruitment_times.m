function recruitment = get_recruitment_times(mea, varargin)
% Compute the recruitment and termination time of each electrode.
% Recruitment time is defined as the first moment with sustained (>= 1 s)
% elevated firing rate (>= 2 sd). Termination time is the first time in the
% last 10 seconds of the seizure in which the firing rate is consistently
% (>= 5s) below the median.
% To return <mea> after computing event_times, provide a variable name
%		e.g. get_recruitment_times(mea, 'mea');


time = mea.Time; time = time();
if ~isfield(mea, 'params'), mea.params = init_mea_params(); end
[~, ~, mea] = mua_events(mea);

if nargin > 1, assignin('caller', varargin{1}, mea); end

fr = zeros(size(mea.mua));
[~, nCh] = size(fr);
fr(mea.event_inds) = 1;  % Store firing times
fr = zscore(smoothdata(fr, 'gaussian', mea.SamplingRate / 2));  % smooth with 500 ms gaussian kernel (merricks, 2016)
peaks = smoothdata(fr > 2, 'movmean', mea.SamplingRate) >= 1;  % Highlight sustained periods of elevated firing
recruitment_time = arrayfun(@(ii) time(find([diff(peaks(:, ii)); 1], 1)), 1:nCh);  % Get time of first peak
recruitment_time(recruitment_time == time(end)) = nan;  % Exclude peaks at the last index
recruitment_time(isoutlier(recruitment_time)) = nan;

EOS_mask = time >= time(end) - mea.Padding(2) - 10;  % Isolate the end of the seizure
tE = time(EOS_mask);
troughs = smoothdata(fr(EOS_mask, :) < quantile(fr, .5), 'movmean', 5*mea.SamplingRate) >= 1;  % highlight sustained silence
termination_time = arrayfun(@(ii) tE(find([diff(troughs(:, ii)); 1], 1)), 1:nCh);  % Get time of first peak
termination_time(termination_time == time(end)) = nan;  % Exclude peaks at the last index
termination_time(isoutlier(termination_time)) = nan;

% Return results
P = mea.Position;                   
P(mea.BadChannels, :) = [];
addy = sub2ind([10 10], P(:, 1), P(:, 2));
recruitment.time = recruitment_time;
recruitment.pos = P;
recruitment.addy = addy;
recruitment.N_rec = sum(isfinite(recruitment_time));
if recruitment.N_rec > 1, recruitment.rate = recruitment.N_rec / range(recruitment.time); else, recruitment.rate = nan; end

recruitment.term_time = termination_time;
recruitment.term_rate = sum(isfinite(termination_time)) / range(termination_time);
