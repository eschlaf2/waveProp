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
if ~isfield(mea, 'event_inds'), [~, ~, mea] = mua_events(mea); end

if nargin > 1, assignin('caller', varargin{1}, mea); end

fr = zeros(size(mea.mua));
[~, nCh] = size(fr);
fr(mea.event_inds) = 1;  % Store firing times
fr = smoothdata(fr, 'gaussian', mea.SamplingRate / 2) * mea.SamplingRate;  % smooth with 500 ms gaussian kernel
[pks, locs, w, p] = deal(cell(nCh, 1));
for ii = 1:nCh, [pks{ii}, locs{ii}, w{ii}, p{ii}] = findpeaks(fr(:, ii), time, 'MinPeakWidth', 1, 'MinPeakHeight', 2*nanstd(fr(:, ii))); end


recruitment_time = get_first(locs);  % get time of first peak
width = get_first(w);
height = get_first(pks);

mask = isoutlier(recruitment_time);
recruitment_time(mask) = nan;  % exclude outliers
width(mask) = nan;
height(mask) = nan;

% Time when each recruited electrode stops firing for 5 seconds
EOS_mask = time >= max(time(end) - mea.Padding(2) - 30, 10);  % Isolate the end of the seizure
tE = time(EOS_mask);
troughs = smoothdata(fr(EOS_mask, :) < quantile(fr, .5), 'movmean', 5*mea.SamplingRate) >= 1;  % highlight sustained silence
termination_time = arrayfun(@(ii) tE(find([diff(troughs(:, ii)); 1], 1)), 1:nCh);  % Get time of first peak
termination_time(termination_time == time(end)) = nan;  % Exclude peaks at the last index
termination_time(isnan(recruitment_time)) = nan;  % Only include recruited electrodes
termination_time(isoutlier(termination_time)) = nan;  % Exclude outliers

% Return results
P = mea.Position;                   
P(mea.BadChannels, :) = [];
addy = sub2ind([10 10], P(:, 1), P(:, 2));
recruitment.time = recruitment_time;
recruitment.pos = P;
recruitment.addy = addy;
recruitment.N_rec = sum(isfinite(recruitment_time));
recruitment.rate = get_rate(P, recruitment_time);
recruitment.width = width;
recruitment.height = height;
recruitment.sd = nanstd(fr);

recruitment.term_time = termination_time;
recruitment.term_rate = get_rate(P, termination_time);

end

function X = get_first(C)
	 idx = ~cellfun(@isempty, C);
	 X = nan(size(C));
	 X(idx) = cellfun(@(c) c(1), C(idx));
end

function rate = get_rate(P, t)
	[~, idx] = max(pdist([P t(:)]));
	p = nchoosek(1:numel(t), 2);
	rate = abs(pdist(P(p(idx, :), :)) ./ diff(t(p(idx, :))));
end
