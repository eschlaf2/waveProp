function sz = plot_summary_stats(stats, h)
% sz = plot_summary_stats(stats, h=figure())
% Plots a horizontal histogram of requested summary stats


if nargin < 1, stats = []; end
if nargin < 2, h = []; end

[stats, h] = parse_inputs_(stats, h);

N = numel(stats);
r = ceil(N / 4);
c = min(4 * r, N);

sz = BVNY.load_seizures;
sz.display_names = compose('P%d', sz.patientAlt);
sz.style = set_style_(sz);

ii = 1;
for ss = stats
    sz = compute_stat_(sz, ss{:}); 
    barh_(sz, ss{:}, subplot(r, c, ii));
    ii = ii + 1;
end

add_legend_(h.Children(1))
rearrange_axes_(h)

end

%% Local functions
function [stats, h] = parse_inputs_(stats, h)
    valid = dir('summary_stats/');
    if isempty(stats), stats = {'duration', 'numevents', 'spatialcoherence'}; end
    if ischar(stats), stats = {stats}; end
    
    if isempty(h), h = figure('units', 'inches', 'position', [0 0 6 2]); end
    stats = cellfun(@(s) validatestring(s, {valid.name}), stats, 'uni', 0);
    if isnumeric(h), h = figure(h); end
    set(h, 'units', 'normalized');
    set(0, 'currentfigure', h); 
end

function sz = compute_stat_(sz, stat)
    for ii = 1:height(sz)
        value = load(fullfile('summary_stats', stat, ...
            [sz.patient{ii} '_' num2str(sz.seizure(ii))]));
        data = value.(stat);
        if numel(data) > 1, data = nanmedian(data); end
        sz.(stat){ii} = data;
    end
    sz.(stat) = cat(1, sz.(stat){:});
end

function S = summarize_(data)
    S = struct;
%     data = cat(1, sz.(ss{:}){:});
%     data = cat(1, data{:});
    X = 1:length(data);
    [bins, edges] = discretize(data, 10);
    S.counts = histcounts(data, edges);
    S.centers = edges(1:end-1) + diff(edges(1:2)) / 2;
    cumcounts = arrayfun(@(ii) find(X(bins == bins(ii)) == ii), X);
%     T.names = fieldnames(S.(ss{:}));
    S.X = cumcounts;
    S.Y = S.centers(bins);
    S.data = data;
end

function barh_(sz, stat, ax)
    T = summarize_(sz.(stat));
    scale = max(T.counts) / 31;
    barh(ax, T.centers, T.counts, ...
        'facealpha', .3, 'edgecolor', 'none');
    hold(ax, 'on')
%     ii = 1;
%     for f = names
    for ii = 1:height(sz)
%         ind = find(strcmpi(T.names, fix(f{:})));
        plot(ax, ii * scale, T.data(ii), ...
            'displayname', sz.display_names{ii}, ...
            sz.style{ii}{:});
%         ii = ii + 1;
    end
    hold(ax, 'off');
    xlim(ax, [0, 1.1*max(T.X)]);
    title(ax, stat)
    ax.XTick = [1 max(T.X)];
    ax.YTick = ax.YTick([1 end]);
    xlabel(ax, 'Counts')
    
end

function add_legend_(ax)
    ln = flipud(findobj(ax, 'type', 'line'));
    names = get(ln, 'displayname');
    [~, inds] = unique(names);
%     mask = [1; diff(G(:))] ~= 0;
    lgd = legend(ln(sort(inds)));
    lgd.Box = 'off';
    lgd.Position([1 2]) = ...
    [0.02, .5-lgd.Position(4)/2];
end

function rearrange_axes_(h)
    ax = findobj(h, 'type', 'axes');
    shift = .98 - sum(ax(1).Position([1 3]));
    for aa = ax'
        aa.Position(1) = aa.Position(1) + shift;
        aa.OuterPosition([2 4]) = [.02 .96];
    end 
end

function style = set_style_(sz)
G = findgroups(sz.patient);
disagree = sz.disagree > 1;
cmap = lines(max(G));
markers = 'o^';

style = arrayfun(@(g, d) ...
    {'color', cmap(g, :) * d, ...
    'linewidth', 1.5, ...
    'marker', markers(ceil(g / 7)), ...
    'markerfacecolor', cmap(g, :)}, ...
    G, ~disagree, 'uni', 0);
 
end


