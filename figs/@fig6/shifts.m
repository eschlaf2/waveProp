function figs = shifts(F)
    % Makes the plots that show shifts and Φ-progressions. Includes methods
    % for stable interval detections and merges.
        
    figs = [ ...  % figures for fig6
        combine_stable_discrete_figs_(F, 'orientation', 'horizontal') ...
        oneS_stable_modes_(F, "MG63", 4, "M", ...
            'create_fig', true, 'center_by', 'iw').Parent ... %% One example: stable modes in MG63 s4
        oneS_stable_modes_(F, "MG49", 36, "M", ...
            'create_fig', true, 'center_by', 'iw').Parent ... %% One example: stable modes in MG49 s36
        hist_shifts_discrete_dirs_(F) ...  %% Histograms
        allP_nshifts_v_ndirs_(F) ...  Nshifts v. Ndirs
        allP_stable_modes_(F, 'center_by', 'iw') ...  Stable modes in all patients
    ];

    make_words_(get_data_(F));
    
end

%% Local funs

function make_words_(data)

    % Setup: Number metrics so that M is first, make a function to get the
    % stats
    data.Gmtc = double(data.metric == "D10") + 1;
    keys = ["patient" "seizure" "Gmtc"];
    
    % Range of shift magnitudes
    % We detect shifts between XXX° and YYY° (*mean* XXM◆XXD, *qq* [YYM, ZZM]◆[YYD, ZZD], M◆D10)
    msg = ['We detect shifts between %0.0f° and %0.0f° ' ...
        '(median %0.0f◆%0.0f°, IQR [%0.0f, %0.0f]◆[%0.0f, %0.0f]°, M◆D10)\n'];
    
    data.shifts_d = rad2deg(data.shifts);
    maxmin = quantile(data.shifts_d, [0 1]);
    stats = groupsummary(data, 'Gmtc', @(x) quantile(x, [.5 .25 .75]), 'shifts_d');
    
    fprintf(msg, maxmin(1), maxmin(2), ...
        stats.fun1_shifts_d([1 2 3 5 4 6]));
    
    
    % Number of shifts per seizure
    % shifts per seizure (*mean* XXM◆XXD, *qq* [YYM, ZZM]◆[YYD, ZZD], M◆D10)
    msg = ['shifts per seizure (median %d◆%d, ' ...
        'range [%d, %d]◆[%d, %d], M◆D10)\n'];    
    Nshifts = groupsummary(data, keys, 'nnz', 'shifts');
        % ... get the number of shifts in each seizure
    stats = groupsummary(Nshifts, 'Gmtc', @(x) quantile(x, [.5 0 1]), 'nnz_shifts');
        % ... and summarize
    fprintf(msg, stats.fun1_nnz_shifts([1 2 3 5 4 6]));
    
    
    % Number of seizures with multiple shifts
    % multiple shifts occur in many seizures (XXX◆YYY/17 seizures, M◆D10)
    msg = 'multiple shifts occur in many seizures (%d◆%d/17 seizures, M◆D10)\n';
    stats = groupsummary(Nshifts, 'Gmtc', @(x) sum(x > 1), 'nnz_shifts');
    fprintf(msg, stats.fun1_nnz_shifts);
    
    
    % Number of unique phases per seizure
    % We find 0 to 4 unique phases per seizure (range [XXX,YYY]◆[XXX,YYY])
    msg = ['We find %d to %d unique phases per seizure ' ...
        '(median %d◆%d, range [%d, %d]◆[%d, %d])\n'];
    Nphases = groupsummary(data, keys, 'max', 'phase_num_adj');
    minmax = quantile(Nphases.max_phase_num_adj, [0 1]);
    stats = groupsummary(Nphases, 'Gmtc', @(x) quantile(x, [.5 0 1]), 'max_phase_num_adj');
    
    fprintf(msg, minmax, stats.fun1_max_phase_num_adj([1 2 3 5 4 6]))
    
    
    % Number of seizures with repeated phases
    % Of these, AAA◆BBB seizures have repeated phases (n = 17 seizures, M◆D10)
    msg = 'Of these, %d◆%d seizures have repeated phases (_n_ = %d seizures, M◆D10)\n';
    has_return = groupsummary(data, keys, 'max', 'isreturn');
    stats = groupsummary(has_return, 'Gmtc', 'sum', 'max_isreturn');
    
    fprintf(msg, stats.sum_max_isreturn, stats.GroupCount(1))
    
    
    % Number of seizures matching each hypothesis
    % From this we conclude that there are YYY/17 seizures consistent with the fixed source hypothesis, XXX/17 seizures consistent with the IW hypothesis, and  YYY/17 seizures consistent with neither existing hypothesis.
    msg = ['From this we conclude that there are %d◆%d/17 ' ...
        'seizures consistent with the fixed source hypothesis, %d◆%d/17 ' ...
        'seizures consistent with the IW hypothesis, and  %d◆%d/17 ' ...
        'seizures consistent with neither existing hypothesis (M◆D10).\n'];
    
    fs_match_ = @(Nshift, Nphase) Nshift == 0 & Nphase == 1;
    iw_match_ = @(Nshift, Nphase) Nshift == 1 & Nphase == 2;
    
    each_sz = join(Nphases, Nshifts, 'keys', keys);
    each_sz.fs_match = fs_match_(each_sz.nnz_shifts, each_sz.max_phase_num_adj);
    each_sz.iw_match = iw_match_(each_sz.nnz_shifts, each_sz.max_phase_num_adj);
    each_sz.no_match = ~(each_sz.fs_match | each_sz.iw_match);
    
    assert(all(each_sz.iw_match + each_sz.fs_match + each_sz.no_match == 1))  % make sure there's no overlap and all seizures are assigned
    sumry = groupsummary(each_sz, 'Gmtc', 'sum', ["iw_match" "fs_match" "no_match"]);
    fprintf(msg, sumry.sum_fs_match, sumry.sum_iw_match, sumry.sum_no_match)
        
    
end

function phase_num = merge_algo_(S, varargin)
    
    PP = parse_inputs_(varargin{:});
    THRESH = PP.thresh_merge_phases;  
    SHOW_MERGES = PP.show_merges;
    
    phase_num = S.phase_num;
    t0 = S.onset_time;
    tF = S.offset_time;
    phi0 = S.phi0;
    phiF = S.phiF;

    % take the closest within thresh?
    while true
        dphi_all = abs(fix_angle(phi0 - phiF'));
        dt_all = t0 - tF';
        merge_score = dphi_all;
        merge_score(dt_all <= 0 | dphi_all > THRESH) = nan;
        merge_score = 1./ merge_score;

        if ~any(merge_score(:)), break; end
        
        % Merge all if all merge scores are > 2/THRESH
        if all(merge_score(isfinite(merge_score)) > 2/THRESH)
            [b_later, b_earlier] = find(merge_score > 0);
            dest = min(b_earlier);
            assert(all(ismember(b_earlier(b_earlier > dest), b_later)))  % If this isn't true, you have a bug
            phase_num(ismember(phase_num, b_later)) = dest;
            if SHOW_MERGES
                mm = sprintf('%d ', unique(b_later));
                fprintf('Merging all: [%s\b] to %d (%s %d %s)\n', ...
                    mm, dest, S.patient(1), S.seizure(1), S.metric(1));
            end
            break
        end

        [b_later, b_earlier] = ...
            find(merge_score == max(merge_score, [], 'all'), 1);

        tF(b_earlier) = tF(b_later);
        phiF(b_earlier) = phiF(b_later);

        [t0(b_later), tF(b_later)] = deal(nan);
        [phi0(b_later), phiF(b_later)] = deal(nan);
        phase_num(phase_num == b_later) = phase_num(b_earlier);
        
        if SHOW_MERGES
            fprintf('Merged %d to %d (%s %d %s)\n', ...
                b_later, b_earlier, S.patient(1), S.seizure(1), S.metric(1))
        end
    
    end
    
    phase_num = findgroups(phase_num);
        
end

function [sumry, G] = summarize_phases_(data)
    % Summarizes information about data based on grouping variables 
    % ["patient" "seizure" "metric", "phase_num"]. Used to merge phases.
    
    [G, sumry] = findgroups(data(:, ["patient" "seizure" "metric", "phase_num"]));

    sumry.isreturn = false(height(sumry), 1);
    sumry.onset_time = splitapply(@(time) time(1), data.test_times, G);
    sumry.offset_time = splitapply(@(time) time(end), data.test_times, G);
    sumry.phi0 = splitapply(@(phi) phi(1), data.ang_mode, G);
    sumry.phiF = splitapply(@(phi) phi(end), data.ang_mode, G);
    sumry.phase_num_adj = sumry.phase_num;
    sumry.row_num = (1:height(sumry))';
    
end

function sumry = merge_phases_(F, data, varargin)
    % summ = adjust_phases_(F, 'thresh_merge_phases'=deg2rad(2))
    %
    % Creates a table with one row for each (seizure, metric, phase
    % number). Merges phases if the drift rate between the phases is less
    % than <thresh_merge_phases> (and total drift < 60°). Phases need not
    % be consecutive to be merged.
    %
    
    if nargin < 2 || isempty(data) 
        data = get_data_(F);
    elseif ~istable(data)
        varargin = [data varargin];
        data = get_data_(F);
    end
    
    PP = parse_inputs_(varargin{:});
    
    sumry = summarize_phases_(data);

    [Gsz, temp] = findgroups(sumry(:, ["patient" "seizure" "metric"])); %#ok<ASGLU>
    
    % test cases
    if 0
        get_ind_ = @(pat, sz, mtc) ...
            find(temp.patient == pat & temp.seizure == sz & temp.metric == mtc); %#ok<UNRCH>
        tc = [ ...
            get_ind_("CUCX3", 6, "M") ... % [1 2 3 2 1 2]
            get_ind_("CUCX5", 3, "M") ...  % [1 1 1 1]
            get_ind_("CUCX5", 6, "M") ...  % [1 1 1 1]
            get_ind_("MG63", 3, "M") ...  % [1 2 2]
            get_ind_("MG63", 4, "M") ...  % [1 2 3]
            get_ind_("MG49", 36, "M") ... % [1 2 3 2 3 2]
            get_ind_("BW09", 2, "D10") ... % [1 2 3 3 4]
            get_ind_( "CUCX3", 2, "D10") ... % [1 2 1]
            get_ind_("CUCX3", 6, "D10") ... % [1 2 3 1]
            get_ind_("CUCX5", 6, "D10") ... % [1 2 3 2 2]
            get_ind_("MG49", 36, "D10") ... % [1 1 2 1]
            ...
            get_ind_("MG49", 36, "M") ... % [1 2 3 2 3 2]
            ];
        for ii = tc
            mask = Gsz == ii;
            S = sumry(mask, :);
            assert(all(S.phase_num == (1:numel(S.phase_num))'));  % Make sure these are 1:N

            phase_num = merge_algo_(S, PP);
            
            fprintf('%s s%d %s: ', S.patient(1), S.seizure(1), S.metric(1));
            fprintf('%d ', phase_num);
            fprintf('\n');

        end


    end

    for ii = 1:max(Gsz)  % for each seizure+metric

        mask = Gsz == ii;
        if sum(mask) < 2, continue; end

        S = sumry(mask, :);

        assert(all(S.phase_num == (1:numel(S.phase_num))'));  % Make sure these are 1:N
        
        phase_num = merge_algo_(S, PP);
        
        sumry.isreturn(S.row_num) = phase_num < (1:numel(phase_num))';
        sumry.phase_num_adj(S.row_num) = phase_num;
    end

end

function ax = oneS_stable_modes_(F, pat, sz, mtc, varargin)
    
    DEF_PAT = "MG63";
    DEF_SZ = 4;
    DEF_MTC = "M";
    
    if nargin < 4 || isempty(mtc), mtc = DEF_MTC; end
    if nargin < 3 || isempty(sz), sz = DEF_SZ; end
    if nargin < 2 || isempty(pat), pat = DEF_PAT; end
    
    if isaxes(varargin{1}), ax = varargin{1}; varargin(1) = []; end
    PP = parse_inputs_(varargin{:});
    
    % Create a figure for each metric 
    if PP.create_fig
        h = figure('position', [0 0 2 1]*2, ...
            'name', sprintf('%s_%d_%smodes%s', pat, sz, mtc, PP.fig_suffix));
        ax = axes(h);
    elseif ~exist('ax', 'var'), ax = gca;
    end
    
            
    keys = ["patient" "seizure" "metric"];
    % don't use get_data_ here... show all seizures, not just those with IW
    data_c = F.CompiledData;
    data_c.time_rel = data_c.time - fillmissing(data_c.iw_center, 'constant', 0);
    data_c.dir_rel = fix_angle(data_c.dir - data_c.iw_angle);
    
    data_s = F.StableDirs;
    data_s.time_iw_rel = data_s.test_times - fillmissing(data_s.iw_center, 'constant', 0);
    data_s.mode_iw_rel = fix_angle(data_s.ang_mode - data_s.iw_angle);
    data_s = join(data_s, merge_phases_(F, data_s), ...
        'keys', [keys "phase_num"]);
    
    % Get the patient and seizure from the argin
    mask_ = @(tbl) tbl(tbl.patient == pat & tbl.seizure == sz & tbl.metric == mtc, :);
    data_c = mask_(data_c);
    data_s = mask_(data_s);
    
    data = renamevars(outerjoin( ...
        groupsummary(data_c, [keys "iw_angle" "iw_center"]), ...
        groupsummary(data_s, keys), ...
        'keys', keys, 'mergekeys', true), ...
        ["GroupCount_left" "GroupCount_right"], ["raw_count" "stable_count"]);
    
    assert(height(data) == 1);  % this should only have one seizure
    
    raw_times = data_c.time_rel;
    raw_dir = data_c.dir;
    stable_times = data_s.time_iw_rel;
    stable_dir = data_s.ang_mode;
    phase_num = data_s.phase_num;
    phases_adj = data_s.phase_num_adj;
    duration_of_phase = data_s.duration_of_phase;
    iw_angle = fillmissing(data.iw_angle, 'constant', 0);
    
    
    raw_dir_adj = raw_dir;
    stable_dir_adj = stable_dir;
    ylim_offset = 0;
    switch PP.center_by
        case 'longest'
            % Center modes so that longest phase is at direction 0
            dirs = stable_dir;
            dur = duration_of_phase;
            mode_longest = dirs(find(dur == max(dur), 1));
%             raw_dir_adj = fix_angle(raw_dir - mode_longest);
%             stable_dir_adj = fix_angle(stable_dir - mode_longest);
%             iw_adj = fix_angle(data.iw_angle - mode_longest);
            ylim_offset = mode_longest;
        case 'iw'
            % Center modes so that iw dir is at direction 0
%             raw_dir_adj = fix_angle(raw_dir - data.iw_angle);
%             stable_dir_adj = fix_angle(stable_dir - data.iw_angle);
%             iw_adj = 0;
            ylim_offset = iw_angle;
        case 'first'
            % Center modes so that the first dir is at direction 0
            if isempty(stable_dir), first_dir = 0;
            else, first_dir = stable_dir(1);
            end
%             raw_dir_adj = fix_angle(raw_dir - first_dir);
%             stable_dir_adj = fix_angle(stable_dir - first_dir);
%             iw_adj = fix_angle(data.iw_angle - first_dir);
            ylim_offset = first_dir;
        case 'none'
            % do nothing
    end

    switch PP.color_modes_by
        case 'phase'
            cc = phase_num;
            CLIM = [1 inf];
        case 'phase_adj'
            cc = phases_adj;
            CLIM = [1 inf];
        case 'time_rel'
            cc = data_s.time_rel;
            CLIM = [0 1];
        otherwise
            error('Missing definition for parameter <color_modes_by>.') 
    end

    % Unwrap dir_mode to draw connected
    dir_mode_u = unwrap(stable_dir_adj) + offsets4unwrap(unwrap(stable_dir_adj));

    % Get a good size for raw points
    uu = ax.Units;
    set(ax, 'units', 'points');
    raw_point_siz = max(min(ax.Position([3 4])/30), 2);
    ax.Units = uu;
    
    % Show directions
    plot(ax, raw_times, rad2deg(raw_dir_adj + 2*pi*[-1 0 1]), 'k.', ...
        'markersize', raw_point_siz, 'displayname', 'raw');
    yline(ax, rad2deg(iw_angle), 'k--');
    if ~isempty(stable_dir_adj)
        hold(ax, 'on');
        splitapply(@(x, y) plot(x, y, 'color', .15*[1 1 1]), ...
            stable_times, rad2deg(dir_mode_u), findgroups(phase_num));
        scatter(reshape(stable_times + [0 0 0], [], 1), ...
            reshape(rad2deg(stable_dir_adj + 2*pi*[-1 0 1]), [], 1), ...
            raw_point_siz.^2 * 3, reshape(cc + [0 0 0], [], 1), 'filled');
        hold(ax, 'off');
        set(ax, 'colormap', metric_cmap_(mtc), 'clim', CLIM);
    end
    title(ax, sprintf('%s %d', pat, sz));
    set(ax, 'ylim', 180*[-1 1] + rad2deg(ylim_offset));
    yticks(-360:60:360);
    xlims = round(([raw_times(1), raw_times(end)] + [-2.5 2.5])/5) * 5;  % to the nearest 5 that shows all the data
    xlim(ax, xlims);
    box(ax, 'off');

end

function hD = allP_shifts_and_phases_(F, varargin)

    PP = parse_inputs_(varargin{:});
    hS = allP_shifts_(F, PP);
    hD0 = allP_discrete_dirs_(F, PP);
    fig_name = hD0.Name;

    axS = findobj(hS, 'type', 'axes');
    axD = findobj(hD0, 'type', 'axes');

    % Convert units to inches to preserve size
    set([axS axD], 'units', 'inches');

    % Move axes to a bigger figure
    hD = figure('position', hD0.Position + [0 0 2 0], 'name', fig_name);
    axD = copyobj(axD, hD);
    axD(1).OuterPosition(1) = hD.Position(3)/2;  % Move to right half

    % Copy shift axes to new figure
    axSN = copyobj(axS, hD);
    for ii = 1:2
        a_ = axSN(ii); a0_ = axD(ii);
        a_.Position(1) = sum(a0_.Position([1 3]));  % align left
        a_.Position(3) = a0_.Position(3)/2; % match half-width

        % apply formatting to shift axis
        set(a_, 'ytick', a0_.YTick, 'yticklabel', [], 'xtick', 0:60:180);
        set(a_.YAxis, 'ticklength', [0 0]);

        % apply formatting to both axes
        set([a_, a0_], 'xticklabelrotation', 45, 'xminorgrid', 'on')
        set(ax, 'xticklabelrotation', 45, 'xminorgrid', 'on')

        % Set minor grid tick locations
        arrayfun(@(aa) set(aa.XAxis, 'MinorTickValues', -180:15:180), [a_, a0_])
        
        % Match linewidth across axes
        ln0 = findobj(a0_, 'type', 'line');
        ln = findobj(a_, 'type', 'line');
        set(ln, 'linewidth', ln0(1).LineWidth);
        

    end
    
    close([hS hD0]);

end

function hD = combine_stable_discrete_figs_(F, varargin)

    PP = parse_inputs_(varargin{:});
    hS = allP_shifts_(F, PP);
    hD0 = allP_discrete_dirs_(F, PP);
    fig_name = hD0.Name;

    axS = findobj(hS, 'type', 'axes');
    axD = findobj(hD0, 'type', 'axes');

    % Convert units to inches to preserve size
    set([axS axD], 'units', 'inches');

    % Move axes to a bigger figure
    hD = figure('position', hD0.Position + [0 0 2 0], 'name', fig_name);
    axD = copyobj(axD, hD);
    axD(1).OuterPosition(1) = hD.Position(3)/2;  % Move to right half

    % Copy shift axes to new figure
    axSN = copyobj(axS, hD);
    for ii = 1:2
        a_ = axSN(ii); a0_ = axD(ii);
        a_.Position(1) = sum(a0_.Position([1 3]));  % align left
        a_.Position(3) = a0_.Position(3)/2; % match half-width

        % apply formatting to shift axis
        set(a_, 'ytick', a0_.YTick, 'yticklabel', [], 'xtick', 0:60:180);
        set(a_.YAxis, 'ticklength', [0 0]);

        % apply formatting to both axes
        set([a_, a0_], 'xticklabelrotation', 45, 'xminorgrid', 'on')
        set([a_, a0_], 'xticklabelrotation', 45, 'xminorgrid', 'on')

        % Set minor grid tick locations
        arrayfun(@(aa) set(aa.XAxis, 'MinorTickValues', -180:15:180), [a_, a0_])
        
        % Match linewidth across axes
        ln0 = findobj(a0_, 'type', 'line');
        ln = findobj(a_, 'type', 'line');
        set(ln, 'linewidth', ln0(1).LineWidth);
        

    end
    
    close([hS hD0]);
    
    % Change to horizontal orientation if necessary
    if strcmpi(PP.orientation, 'horizontal')
        
        ax = findobj(hD, 'type', 'axes');
        % Rotate all the labels
        set(ax, ...
            'positionconstraint', 'innerposition', ...
            'yticklabelrotation', -45, ...
            'xticklabelrotation', -90, ...
            'xaxislocation', 'top');
        for aa = ax', aa.Position(2) = .1; end

        % correct offsets to advance positive
        correct_yy_ =@(ln) set(ln, 'YData', 2*ceil(ln.YData) - ln.YData);

        ln = [findobj(ax, 'type', 'line'); ...
            findobj(ax, 'type', 'scatter')];
        arrayfun(@(ll) correct_yy_(ll), ln);

    end

end

function gap_counting_details_(F, varargin)

    PP = parse_inputs_(varargin{:});

    % Parameters
    THRESH1 = PP.thresh_discrete;
    THRESH2 = PP.thresh_shift;

    % Definitions
    W = 3;  % One sample on either side. 
    isoutlier_ = @(x) ...
        ( (x - movmedian(x, W)) ./ movmad(x, W) ) > 3 ...
        & movmad(x, W) > 0 ...
        & x <= pi;  % this isn't necessary, but I want to see non-wraparound outliers
    isgap_ = @(diffs) ...
        diffs > THRESH1 ...
        | ( isoutlier_(diffs) & diffs > THRESH2 );


    % Load data
    data = get_data_(F, PP);


    % Compute gaps for each seizure and metric
    [G, pat, sz, mtc] = findgroups(data.patient, data.seizure, data.metric);
    Ngaps = nan(max(G), 1);
    Nmodes = Ngaps;
    [diffs, outs, md] = deal(cell(numel(Ngaps), 1));
    d3 = nan(numel(Ngaps), 3);
    for iG = unique(G)'

        dirs = data.ang_mode(G == iG);  % get all modes
        dirsS = sort(dirs);  % ... and sort them

        % Compute differences between consecutive angles
        diffs{iG} = abs(diff([dirsS; dirsS(1) + 2*pi]));

        % Extra details for verification/inspection
        % ... The three largest shifts for each group;
        dd = sort(diffs{iG}(diffs{iG} <= pi), 'desc');
        Nd = min(numel(dd), 3);
        d3(iG, 1:Nd) = dd(1:Nd);

        % ... size of outliers and local mean
        outs{iG} = diffs{iG}(isoutlier_(diffs{iG}));
        temp_ = movmedian(diffs{iG}, W);
        md{iG} = temp_(isoutlier_(diffs{iG}));

        % Count the number of gaps
        num_gaps = sum( isgap_(diffs{iG}) );  % Count the number of gaps
        Ngaps(iG) = num_gaps;

        % ... and use it compute the number of sources: Nmodes = Ngaps except
        % possibly in the zero case (modes wrap around the full circle without
        % a gap). I think it's best here to leave this value at zero so that we
        % are counting *discrete source directions*. I also don't think this
        % comes up in this dataset anyway...
        Nmodes(iG) = num_gaps;  

    end

    % Compile all diffs into a vector
    diffs_mat = cat(1, diffs{:});

    % Show details
    hDetail = figure;
    T = BVNY.tiledlayout(hDetail, 'flow');

    % ... hist of all diffs
    ax = nexttile(T);
    edges = [0:15 30:15:180];
    histogram(ax, rad2deg(diffs_mat), [edges inf])
    xticks(ax, 0:15:180);
    xlim(ax, [0 195])
    ylim(ax, [0 100])
    title('Diffs between modes')

    % ... three largest diffs
    ax = nexttile(T);
    lnM = plot(ax, rad2deg(d3(mtc == "M", :)'), 'k-o', 'displayname', 'M');
    hold(ax, 'on');
    lnD10 = plot(ax, rad2deg(d3(mtc == "D10", :)'), 'r-o', 'displayname', 'D10');
    hold(ax, 'off');
    ylabel(ax, 'Diffs [°]')
    title(ax, 'Three largest diffs');
    legend(ax, [lnM(1), lnD10(1)]);

    % ... size of outliers vs. surrounding mean
    ax = nexttile(T);
    xo = rad2deg(cat(1, md{:}));
    yo = rad2deg(cat(1, outs{:}));
    jitter_ =@() (rand(size(xo)) -.5) * .5;
    plot(ax, xo + 2*jitter_(), yo + jitter_(), '.');
    QQ = quantile([xo(:); yo(:)], [0 1]);
    line(ax, QQ, QQ, 'color', 'k', 'linestyle', '--')
    xlabel(ax, 'local mean diff [°]')
    ylim(ax, [0 rad2deg(THRESH2)] + [-1 1])
    ylabel(ax, 'size of outlier [°]')
    title(ax, sprintf('Outliers lower than THRESH2 (N=%d)', ...
        sum(yo < rad2deg(THRESH2))))

    % ... size of outliers vs. max gap in group
    ax = nexttile(T);
    maxGap = rad2deg(d3(:, 1));
    maxG_ind = cell2mat(arrayfun(@(ii) repmat(ii, 1, numel(outs{ii})), 1:size(d3, 1), 'uni', 0));
    xo = maxGap(maxG_ind);
    yo = rad2deg(cat(1, outs{:}));
    jitter_ =@() (rand(size(xo)) -.5) * .5;
    plot(ax, xo + 2*jitter_(), yo + jitter_(), '.');
    QQ_ = @(x) quantile(x, [0 1]);
    line(ax, QQ_(xo), QQ_(yo), 'color', 'k', 'linestyle', '--')
    xlabel(ax, 'Max gap for group [°]')
    ylim(ax, [0 rad2deg(THRESH2)] + [-1 1])
    ylabel(ax, 'size of outlier [°]')
    title(ax, sprintf('Outliers lower than THRESH2 (N=%d)', ...
        sum(yo < rad2deg(THRESH2))))


    % Show ladder plot of number of gaps
    h = figure('name', sprintf('gaps%s', PP.fig_suffix), ...
        'position', [0 0 3 4]);
    ax = axes(h);
    xx = (mtc == F.Metrics(2)) + 1;
    yy = Ngaps;
    offset = rescale(findgroups(pat, sz), 0, .5) - .25;
    F.gscatter_pat(ax, xx + offset/4, yy + offset, pat, 'seizure', sz);
    legend('location', 'eastout')
    xticks([1 2])
    xticklabels(F.Metrics);
    xlim([0.5 2.5])
    yticks(1:max(yy))
    ylabel('Gaps between modes [#]')




    ax = gca;
    cla

    hold(ax, 'on')
    arrayfun(@(mm) histogram(Nmodes(mtc == mm), 'displayname', mm, ...
            'facecolor', F.Style.(mm).color, 'edgecolor', .85*[1 1 1]), ...
        F.Metrics)
    hold(ax, 'off');
    legend(ax, 'location', 'eastout')
    xlabel('Discrete dirs per seizure [#]');
    ylabel('Counts');

end

function h = hist_shifts_discrete_dirs_(F, varargin)
    
    PP = parse_inputs_(varargin{:});

    % Get shifts and summary data
    data = get_data_(F, PP);
    datS = data(data.isshift, :);
    datJ = summarize_seizures_(F, PP);

    % Create the figure
    h = figure('name', sprintf('hist_shifts_discrete_dirs%s', PP.fig_suffix), ...
        'position', [0 0 6 2]);
    T = F.tiledlayout(h, 1, 3);

    % Convenience functions
    set_bar_ =@(bar) arrayfun(@(bb, mm) ...
        set(bb, 'facecolor', F.Style.(mm).color, 'displayname', mm, ...
        'barwidth', 1.5, 'facealpha', .8, 'edgecolor', 'none'), bar, F.Metrics);
    edge2ctr_ = @(edges) mean([edges(1:end - 1); edges(2:end)]);
    Gmtc_ = @(tbl) (tbl.metric == F.Metrics(2)) + 1;


    % Histogram of all shift magnitudes
    ax = nexttile(T, 1);
        % ... get the data
    Gm = Gmtc_(datS);
    edges = 0:15:180;
    xi = edge2ctr_(edges);
        % ... and make the bar plot
    counts = splitapply(@(x) histcounts(x, edges), rad2deg(datS.dphi), Gm);
    b = bar(ax, xi, counts');
    set_bar_(b);
        % ... and prettify
    xlabel('Magnitude [°]')
    xticks(0:30:180)
    xlim([0 180])
    ylabel('Shifts [#]')


    % Histogram of number of shifts per seizure
    ax = nexttile(T, 2);
        % ... get the data
    Gm = Gmtc_(datJ);
    shifts = fillmissing(datJ.Nbig, 'constant', 0);
    edges = -.5 : max(shifts) + 1;
    counts = splitapply(@(N) histcounts(N, edges), shifts, Gm);
        % ... plot and prettify
    b = bar(ax, edge2ctr_(edges), counts');
    set_bar_(b);
    xlabel('Shifts [#]');
    ylabel('Seizures [#]');


    % Histogram of discrete phases
    ax = nexttile(T, 3);
    Gm = Gmtc_(datJ);
    NDD = datJ.Nphases;
    edges = (min(NDD):max(NDD) + 1) - .5;
    counts = splitapply(@(N) histcounts(N, edges), NDD, Gm);
    b = bar(ax, counts');
    set_bar_(b);
    xlabel('Discrete phases [#]');
    ylabel('Seizures [#]');

end

function h = allP_nshifts_v_ndirs_(F, varargin)
    
    
    PP = parse_inputs_(varargin{:});
    data = summarize_seizures_(F, PP);
    data.Nsmall(data.Nsmall == 0) = nan;

    h = figure('name', sprintf('nshifts_v_ndirs%s', PP.fig_suffix), ...
        'position', [0 0 6 5]);
        
            
    % Get the plotting values and store in a table
    plt = sortrows( ...
        renamevars( stack(data, ["Nbig" "Nsmall"]), ...
            ["Nbig_Nsmall_Indicator", "Nbig_Nsmall", "Nphases"], ...
            ["issmall", "xx", "yy"]), ...
        ["issmall", "metric"], {'desc', 'desc'});
    
    ax_nvn = axes(h);
    plot_nvn_(F, plt, ax_nvn);
    lgd = ax_nvn.Legend;
    set(lgd, 'units', 'inches');
    lgd.Position(1:2) = [0.1 h.Position(4) - lgd.Position(4) - .1]; 
    set(ax_nvn, 'units', 'inches', 'DataAspectRatioMode', 'auto', ...
        'PlotBoxAspectRatioMode', 'auto'); 
    R = 4/ax_nvn.OuterPosition(3);
    ax_nvn.OuterPosition([3 4]) = ax_nvn.OuterPosition([3 4]) * R;
    ly = range(ax_nvn.YLim); lx = range(ax_nvn.XLim);
    ax_nvn.InnerPosition(4) = ly/lx * ax_nvn.InnerPosition(3);
        
    POS = ax_nvn.InnerPosition;
    
    % Create the histograms
    fig_temp = hist_shifts_discrete_dirs_(F, PP);
    ax_src = findobj(fig_temp, 'type', 'axes');
    
    % Move shifts hist into a tile here
    axS = copyobj(ax_src(2), h);
    set(axS, 'units', 'inches', ...
        'innerposition', [POS(1), sum(POS([2 4])), POS(3), POS(4)/2]);
    linkaxes([ax_nvn axS], 'x');
    ylim(axS, [0 inf]);
    
    
    % Show discrete dirs hist
    axD = copyobj(ax_src(1), h);
    set(axD.Children, 'horizontal', 'on');
    xlim(axD, [0 inf]);
    ylim(axD, ax_nvn.YLim);
    set(axD, 'units', 'inches', ...
        'innerposition', [sum(POS([1 3])), POS(2), POS(3)/2, POS(4)]);
    
    
%     axD.InnerPosition([2 4]) = ax_nvn.InnerPosition([2 4]);
    
    close(fig_temp);
    
    % Add numeric labels to vertical bars
    xtips_ = @(b) b.XEndPoints;
    ytips_ = @(b) b.YEndPoints;
    labels_ = @(b) string(b.YData);
    label_bars_ = @(b) text(axS, xtips_(b), ytips_(b), labels_(b), ...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    arrayfun(@(b) label_bars_(b), axS.Children);
    
    % Add numeric labels to horizontal bars
    xtips_ = @(b) b.YEndPoints+.05;
    ytips_ = @(b) b.XEndPoints;
    labels_ = @(b) string(b.YData);
    label_bars_ = @(b) text(axD, xtips_(b), ytips_(b), labels_(b), ...
        'HorizontalAlignment','left',...
        'VerticalAlignment','middle');
    arrayfun(@(b) label_bars_(b), axD.Children);
    
    set([axS axD], 'visible', 'off');
    
end

function plot_nvn_(F, plt, ax)
    
    % Functions to spread out points with same x- and y-values
    G_ = @(x, y) findgroups(x, y);  % find points where x and y overlap
    number_ = @(G) arrayfun(@(ii) ...
        sum(G(1:ii-1) == G(ii)), (1:numel(G))');
    rescale_ = @(x) x * .1;
    spread_ = @(x, y) rescale_(number_(G_(x, y)));
    
    % ... add spread
    spread = rescale(spread_(plt.xx, plt.yy)) * .6;
    for ff = ["xx" "yy"], plt.(ff) = plt.(ff) + spread; end
    
    
    % Draw lines to indicate where each shift goes to a distinct
    % direction
    LL = [-2, max(plt.yy) + 2];
    line(ax, LL, LL + 1.7, 'linestyle', '--', 'color', .15*[1 1 1]);
    line(ax, LL, LL + .3, 'linestyle', '--', 'color', .15*[1 1 1]);
                
    % Make the plots for each metric
    hold(ax, 'on')
    for mm = unique(plt.metric)' 
                
        % Filter the data to metric mm
        plt_ = plt(plt.metric == mm & plt.issmall == "Nbig", :);  % only show large shifts
        
        ln = F.gscatter_pat(ax, plt_.xx, plt_.yy, plt_.patient);
        
%         % Plot *small* and *large* shifts separately
%         G = findgroups(plt_.issmall == "Nbig");  % this makes *small* the first group
%         ln = splitapply(@(xx, yy, pat) ...
%                 F.gscatter_pat(ax, xx, yy, pat)', ...
%             plt_(:, ["xx" "yy" "patient"]), G);
%         
% 
%         % Fade color in *small* shifts
%         fade_ =@(ll) mean([ll.Color; 1 1 1]);
%         for ll = ln(1, :)  
%             set(ll, 'markerfacecolor', fade_(ll), 'color', fade_(ll));
%         end
%          
        
        % Style by metric
        if mm == "D10"
            set(ln, 'markerfacecolor', [1 1 1]);
        else
            set(ln, 'color', [1 1 1], 'linewidth', .5, ...  % thin white border
                'markersize', ln(1).MarkerSize + 1);  %  ... size compensated
        end
        
        % Set Zdata so marker overlap looks good
        arrayfun(@(ll) set(ll, 'zdata', -ll.YData, 'linewidth', 1), ln);
                
    end
    hold(ax, 'off')
    
    % Prettify
    axis(ax, 'equal');
    xlabel('Shifts [#]')
    xlim(ax, [-.5 max(plt.xx) + .5]); 
    xticks(ax, 0:max(plt.xx));
    ylim(ax, quantile(plt.yy, [0 1]) + [-1 1]*.5); 
    yticks(ax, 1:max(plt.yy));
    ylabel('Unique phases [#]')
    legend(ax, 'location', 'northout', 'numcolumns', 3)
    
    % switch x and y axes
%     view(ax, 90, -90);
%     arrayfun(@(ll) set(ll, 'zdata', -ll.ZData), findobj(ax, 'type', 'line'));

end

function data = get_data_(F, varargin)
    % Get data and mask it to stable phases and seizures with IW
    
    PP = parse_inputs_(varargin{:});
    keys = ["patient", "seizure", "metric"];
    data = F.StableDirs;
    data.isshift = logical(data.isshift);
    
    data = join(data, merge_phases_(F, data, PP), ...
        'keys', [keys, "phase_num"]);
    data.isshift = data.isshift & [false; diff(data.phase_num_adj) ~= 0];
    
    % Collapse the phase numbers
    temp = grouptransform(data, keys, ...
        @(x) cumsum(x) + 1, 'isshift');  
    data.phase_num = temp.isshift;
    
    % ... and fix <duration_of_phase>
    temp = grouptransform(data, keys, @range, 'test_times');
    data.duration_of_phase = temp.test_times;
    
    % Mask the data and get yy-values
    mask_ = data.nchannels >= F.MinIWElectrodes;
    data = data(mask_, :);
    data.patientAlt = arrayfun(@(pat) F.PatientNum(pat), data.patient);
    data.yy = findgroups(data.patientAlt, data.seizure);  % keep yy stable for both TOA methods
    
    % shifts was dphi; now dphi/dt
    data.shifts = abs(fix_angle([nan; diff(data.ang_mode)]));
    data.shifts(~data.isshift) = nan;

end

function sumry = summarize_seizures_(F, varargin)
    
    PP = parse_inputs_(varargin{:});
    THRESH = PP.thresh_small_shifts; 
    keys = ["patient", "seizure", "metric"];
    
    
    % Get shifts
    data = get_data_(F, PP);
    datS = data;
    datS.lt_thresh = datS.shifts < THRESH;
    

    % Join data to get counts for the same patients and seizures
    [G, datS_summ] = findgroups(datS(:, keys));
    datS_summ.Nphases = splitapply(@max, datS.phase_num_adj, G);
    datS_summ.Nshifts = splitapply(@sum, datS.isshift, G);
    datS_summ.Nshifts_lt_thresh = splitapply(@sum, datS.lt_thresh, G);
    sumry = datS_summ;
    sumry.Nshifts = fillmissing(sumry.Nshifts, 'constant', 0);
    
    
    for ff = ["Nshifts", "Nshifts_lt_thresh"]
        sumry.(ff) = fillmissing(sumry.(ff), 'constant', 0); 
    end
    
    sumry.Nsmall = sumry.Nshifts_lt_thresh;
    sumry.Nbig = sumry.Nshifts - sumry.Nsmall;
    

end

function cmap = metric_cmap_(mm)
    color = BVNY.Style.(mm).color;
    cmap = make_diverging_colormap( ...
        [color; .15*[1 1 1]], ...
        .7*[1 1 1]);
end

function h = allP_stable_modes_(F, varargin)
    
    PP = parse_inputs_(varargin{:});
            
    data_c = F.CompiledData;  % Just loading this to get all seizures
    [G, pat, sz] = findgroups(data_c.patient, data_c.seizure);
    Gpat = findgroups(pat);
    sznum = cell2mat(splitapply(@(sz) {(1:numel(sz))'}, sz, Gpat));

    rows = max(Gpat);
    cols = max(histcounts(categorical(pat)));

    tile_num_ = @(Gpat, sznum) sub2ind([cols, rows], sznum, Gpat);

    h = gobjects(1, numel(F.Metrics));

    for iM = 1:numel(F.Metrics)

        metric = F.Metrics(iM);  % Look at each metric separately
        
        % Create a figure for each metric 
        h(iM) = figure('position', [0 0 cols*2 rows*1], ...
            'name', sprintf('%sexamples%s', metric, PP.fig_suffix));
        T = BVNY.tiledlayout(h(iM), rows, cols);
        title(T, metric);


        % Plot each seizure in a tile
        for ii = 1:max(G)

            % Get a tile and plot
            ax = nexttile(T, tile_num_(Gpat(ii), sznum(ii)));
            oneS_stable_modes_(F, pat(ii), sz(ii), metric, ax, PP);
        end
    end

end

function h = allP_discrete_dirs_horiz_(F, varargin)
    
    % Parameters
    PP = parse_inputs_(varargin{:});
    ORIENTATION = 'horizontal';
    fig_sfx = PP.fig_suffix;
    
    switch PP.cmap 
        case 'metric'
            cmap_ =@(mm) metric_cmap_(mm);
        case 'red'
            cmap_ = @(mm) make_diverging_colormap( ...
                [[1 0 0]; [0 0 0]], .7*[1 0 0]);
        otherwise
            error('CMAP ''%s'' not recognized, CMAP');
    end
    data = get_data_(F, PP);
    data = join(data, groupsummary( ...
                data(:, ["yy" "phase_num" "phase_num_adj"]), "yy", ...
                'max', ["phase_num" "phase_num_adj"]));

    % Limit which modes to show
    switch PP.modes2show
        case 'onsets'
            % THE FOLLOWING LIMITS THE DATA TO ONLY THE ONSETS OF THE STABLE PHASES.
            % UNCOMMENT TO USE THIS. WITHOUT IT, ALL DETECTED MODES ARE SHOWN
            % Get row numbers of stable onsets 
            G = findgroups(data.patient, data.seizure, data.metric, data.phase_num);
            rows = splitapply(@(row_num) row_num(1), (1:height(data))', G);

            % Restrict data to rows with stable onsets
            data = data(rows, :);
    end

    % Center the directions
    switch PP.center_by
        case 'longest'
            % Center modes so that longest phase is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            mode_longest = splitapply(@(mode, dur) mode(find(dur == max(dur), 1)), ...
                data.ang_mode, data.duration_of_phase, G);
            data.mode_adj = fix_angle(data.ang_mode - mode_longest(G));
            data.iw_adj = fix_angle(data.iw_angle - mode_longest(G));
        case 'iw'
            % Center modes so that iw dir is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            iw_dir = splitapply(@(iw) iw(1), ...
                data.iw_angle, G);
            data.mode_adj = fix_angle(data.ang_mode - iw_dir(G));
            data.iw_adj = 0*data.iw_angle;
        case 'first'
            % Center modes so that the first dir is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            first_dir = splitapply(@(mode) mode(1), ...
                data.ang_mode, G);
            data.mode_adj = fix_angle(data.ang_mode - first_dir(G));
            data.iw_adj = fix_angle(data.iw_angle - first_dir(G));
        case 'none'
            data.mode_adj = data.ang_mode;
            data.iw_adj = data.iw_angle;
    end
    
    switch PP.color_modes_by
        case 'phase'
            % Rescale phase_num to use with cmap
            cc_ =@(temp_dat) (temp_dat.phase_num - 1) ./ ...
                max(temp_dat.max_phase_num - 1, 1);
        case 'phase_adj'
            % Rescale phase_num to use with cmap
            cc_ =@(temp_dat) (temp_dat.phase_num_adj - 1) ./ ...
                max(temp_dat.max_phase_num_adj - 1, 1);
        case 'time_rel'
            cc_ =@(temp_dat) temp_dat.time_rel;
        otherwise
            error('Missing definition for parameter <color_modes_by>.')     
    end

    switch ORIENTATION
        case 'horizontal'
            pos = [0 0 6 4];
            subplot_ = @(ii) subplot(2, 1, ii);
            
            % connecting lines
            connecting_lines_ =@(ax, angu, tt, offset, varargin) ...
                splitapply(@(xx, yy) { ...
                    plot(ax, xx, yy, 'zdata', xx - 100, varargin{:}) ...
                }, tt + offset, angu, findgroups(tt)); 
            
            mode_markers_ = @(ax, dir, tt, offset, cc, Gdot) ...
                splitapply(@(x, y, cc) ...
                    scatter(ax, x, y, 9, cc, 'filled', 'zdata', y), ...
                tt + offset, dir, cc, Gdot);
            
            iw_indicator_ =@(ax, iw_ang, tt, varargin) ...
                splitapply(@(xx, yy) { ...
                    plot(ax, xx', yy', varargin{:}) ...
                }, tt + [0 .6], iw_ang + [0 0], findgroups(tt));
            
            prettify_ = @(ax, labels, tt) { ...
                set(ax, 'ylim', [-1 1] * 200, ...
                'ytick', (-180:60:180), ...
                'xtick', unique(tt), ...
                'xlim', [0 max(tt)+1], ...
                'xticklabel', labels), ...
                xlabel(ax, '[\circ]'), ...
                legend(ax, 'off') ...
                };
            
        case 'vertical'
            pos = [0 0 6 4];
            subplot_ = @(ii) subplot(1, 2, ii);
            
            % connecting lines
            connecting_lines_ =@(ax, angu, tt, offset, varargin) ...
                splitapply(@(xx, yy) { ...
                plot(ax, xx, yy, 'zdata', yy - 100, ...
                varargin{:}) ...
                }, angu, tt - offset, findgroups(tt)); 
            
            % ... add mode markers
            mode_markers_ = @(ax, dir, tt, offset, cc, Gdot) ...
                splitapply(@(x, y, cc) ...
                    scatter(ax, x, y, 9, cc, 'filled', 'zdata', y), ...
                dir, tt - offset, cc, Gdot);
            
            iw_indicator_ =@(ax, iw_ang, tt, varargin) ...
                splitapply(@(xx, yy) { ...
                    plot(ax, xx', yy', varargin{:}) ...
                }, iw_ang + [0 0], tt - [0 .6], findgroups(tt));
            
            
            
            prettify_ = @(ax, labels, tt) { ...
                set(ax, 'xlim', [-1 1] * 200, ...
                'xtick', (-180:60:180), ...
                'ytick', unique(tt), ...
                'ylim', [0 max(tt)+1], ...
                'yticklabel', labels, ...
                'ydir', 'reverse') ...
                xlabel(ax, '[\circ]'), ...
                legend(ax, 'off') ...
                };
            
    end
    
            
    % Show modes for each patient, seizure, and metric
    h = figure('name', sprintf('stable%s', fig_sfx), ...
        'position', pos);
    for ii = 1:2  
        mm = F.Metrics(ii);

        % get the data
        temp_dat = data(data.metric == mm, :);
        
        % Get the plotting values
        dir = rad2deg(temp_dat.mode_adj);
        tt = temp_dat.yy;
        pat = temp_dat.patient;
        sz = temp_dat.seizure;
        iw_ang = rad2deg(temp_dat.iw_adj);

        % Unwrap xx to draw lines connecting phases
        angu = unwrap(temp_dat.mode_adj);  
        angu = rad2deg(angu + offsets4unwrap(angu));

        % offset the x/y-values according to the phase number
        offset = rescale(temp_dat.phase_num) * .6;  % discrete phases
%         offset = temp_dat.time_rel/2;  % smooth flow

        % Color dots by time relative to first and last detected discharges
        [Gdot, ~] = findgroups(temp_dat(:, ["yy", "phase_num"]));
        

        %%
        % plot
        ax = subplot_(ii);
        hold(ax, 'on');
        ln_line = connecting_lines_(ax, angu, tt, offset, ...
            'linewidth', .5, 'color', .15*[1 1 1]); %#ok<NASGU>

        % ... add mode markers
        mode_markers_(ax, dir, tt, offset, cc_(temp_dat), Gdot);

        % ... add IW direction indicator
        iw_indicator_(ax, iw_ang, tt, 'linewidth', .5, 'color', .15*[1 1 1]);
        hold(ax, 'off');
        set(ax, 'colormap', cmap_(mm), 'clim', [0 1]);
        
        % Prettify
        title(sprintf('Stable directions (%s)', mm))
        labels = splitapply(@(pat, sz) sprintf("%s %d", pat(1), sz(1)), ...
            pat, sz, findgroups(tt));
        prettify_(ax, labels, tt);

    end

end

function h = allP_discrete_dirs_(F, varargin)
    
    % Parameters
    PP = parse_inputs_(varargin{:});
    fig_sfx = PP.fig_suffix;
    
    switch PP.cmap 
        case 'metric'
            cmap_ =@(mm) metric_cmap_(mm);
        case 'red'
            cmap_ = @(mm) make_diverging_colormap( ...
                [[1 0 0]; [0 0 0]], .7*[1 0 0]);
        otherwise
            error('CMAP ''%s'' not recognized, CMAP');
    end
    data = get_data_(F, PP);

    % Limit which modes to show
    switch PP.modes2show
        case 'onsets'
            % THE FOLLOWING LIMITS THE DATA TO ONLY THE ONSETS OF THE STABLE PHASES.
            % UNCOMMENT TO USE THIS. WITHOUT IT, ALL DETECTED MODES ARE SHOWN
            % Get row numbers of stable onsets 
            G = findgroups(data.patient, data.seizure, data.metric, data.phase_num);
            rows = splitapply(@(row_num) row_num(1), (1:height(data))', G);

            % Restrict data to rows with stable onsets
            data = data(rows, :);
    end

    % Center the directions
    switch PP.center_by
        case 'longest'
            % Center modes so that longest phase is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            mode_longest = splitapply(@(mode, dur) mode(find(dur == max(dur), 1)), ...
                data.ang_mode, data.duration_of_phase, G);
            data.mode_adj = fix_angle(data.ang_mode - mode_longest(G));
            data.iw_adj = fix_angle(data.iw_angle - mode_longest(G));
        case 'iw'
            % Center modes so that iw dir is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            iw_dir = splitapply(@(iw) iw(1), ...
                data.iw_angle, G);
            data.mode_adj = fix_angle(data.ang_mode - iw_dir(G));
            data.iw_adj = 0*data.iw_angle;
        case 'first'
            % Center modes so that the first dir is at direction 0
            G = findgroups(data.patient, data.seizure, data.metric);
            first_dir = splitapply(@(mode) mode(1), ...
                data.ang_mode, G);
            data.mode_adj = fix_angle(data.ang_mode - first_dir(G));
            data.iw_adj = fix_angle(data.iw_angle - first_dir(G));
        case 'none'
            data.mode_adj = data.ang_mode;
            data.iw_adj = data.iw_angle;
    end


    % Show modes for each patient, seizure, and metric
    h = figure('name', sprintf('stable%s', fig_sfx), ...
        'position', [0 0 4 6]);
    for ii = 1:2  
        mm = F.Metrics(ii);

        % get the data
        temp_dat = data(data.metric == mm, :);
        
        % Get the plotting values
        xx = rad2deg(temp_dat.mode_adj);
        yy = temp_dat.yy;
        pat = temp_dat.patient;
        sz = temp_dat.seizure;
        iw_ang = rad2deg(temp_dat.iw_adj);

        % Unwrap xx to draw lines connecting phases
        xu = unwrap(temp_dat.mode_adj);  
        xu = rad2deg(xu + offsets4unwrap(xu));

        % offset the y-values according to the phase number
        offset = rescale(temp_dat.phase_num) * .6;  % discrete phases
%         offset = temp_dat.time_rel/2;  % smooth flow

        % Color dots by time relative to first and last detected discharges
        [Gdot, ~] = findgroups(temp_dat(:, ["yy", "phase_num"]));
        
        switch PP.color_modes_by
            case 'phase'
                % Rescale phase_num to use with cmap
                temp_dat = join(temp_dat, groupsummary( ...
                    temp_dat(:, ["yy" "phase_num"]), "yy", 'max', "phase_num"));
                cc = (temp_dat.phase_num - 1) ./ ...
                    max(temp_dat.max_phase_num - 1, 1);

            case 'phase_adj'
                % Rescale phase_num to use with cmap
                temp_dat = join(temp_dat, groupsummary( ...
                    temp_dat(:, ["yy" "phase_num_adj"]), "yy", 'max', "phase_num_adj"));
                cc = (temp_dat.phase_num_adj - 1) ./ ...
                    max(temp_dat.max_phase_num_adj - 1, 1);
        
            case 'time_rel'
                cc = temp_dat.time_rel;
            otherwise
                error('Missing definition for parameter <color_modes_by>.')     
        end

        % plot
        ax = subplot(1, 2, ii);
        hold(ax, 'on');
        % ... add connecting line
        ln_line = splitapply(@(xx, yy) { ...
            plot(ax, xx, yy, 'zdata', yy - 100, ...
            'linewidth', .5, 'color', .15*[1 1 1]) ...
            }, xu, yy - offset, findgroups(yy)); %#ok<NASGU>
        % ... add mode markers
        splitapply(@(x, y, cc) ...
            scatter(ax, x, y, 9, cc, 'filled', 'zdata', y), ...
            xx, yy - offset, cc, Gdot);
        % ... add IW direction indicator
        ln_iw = splitapply(@(xx, yy) { ...
            plot(ax, xx', yy', ... 'zdata', yy + 100, ...
            'linewidth', .5, 'color', .15*[1 1 1]) ...
            }, iw_ang + [0 0], yy - [0 .6], findgroups(yy)); %#ok<NASGU>
        
        hold(ax, 'off');
        set(ax, 'colormap', cmap_(mm), 'clim', [0 1]);
        

        % Prettify
        title(sprintf('Stable directions (%s)', mm))
        xlim([-1 1] * 200)
        xticks(-180:60:180)
        xlabel('[\circ]')
        yticks(unique(yy))
        ylim([0 max(yy)+1])
        lbls = splitapply(@(pat, sz) sprintf("%s %d", pat(1), sz(1)), ...
            pat, sz, findgroups(yy));
        yticklabels(lbls)
        legend('off');
        axis ij
    end


    ax = findobj(h, 'type', 'axes');
    linkaxes(ax, 'y');
    ylim(ax, [0 max(data.yy) + 1])
end

function h = allP_shifts_(F, varargin)
    
    % Load the data and get the shifts
    PP = parse_inputs_(varargin{:});
    data = get_data_(F, PP);
    dat_shifts = data(data.isshift, :);   
    
    % Make the figures
    h = figure('name', sprintf('shifts%s', PP.fig_suffix), ...
        'position', [0 0 4 6]);
    
    for ii = 1:2
        mm = F.Metrics(ii);
        temp_dat = dat_shifts(dat_shifts.metric == mm, :);

        figure(h);
        ax = subplot(1, 2, ii);
        xx = rad2deg(abs(temp_dat.dphi));
        yy = temp_dat.yy;
        offset = rescale(temp_dat.phase_num) * .6;  % discrete phases
%         yy = findgroups(temp_dat.patnum, temp_dat.sz);
        pat = temp_dat.patient;
        sz = temp_dat.seizure;
        plot(ax, xx' .* [-1; 1], (yy - offset)'.*[1; 1], 'k-', 'tag', sprintf('%sshift_ln', mm));
        hold(ax, 'on')
        plot(ax, xx, yy - offset, 'o', ...
            'color', F.Style.(mm).color, 'markersize', 3, ...
            'markerfacecolor', F.Style.(mm).color, 'tag', sprintf('%sshift_dot', mm));
        hold(ax, 'off');
        xticks(0:30:180)
        xlim([0 180])
        yticks(unique(yy))
        lbls = arrayfun(@(pat, sz) sprintf("%s %d", pat, sz), pat, sz);
        lbls = arrayfun(@(yu) unique(lbls(yy == yu)), unique(yy));
        yticklabels(lbls);
        legend('off')
        ylim([0, numel(lbls)+1])
        title(sprintf('Shifts (%s)', F.Metrics(ii)))
        xlabel('Magnitude (°)');
        axis ij

        

    end

    % Link y-axes in plots of shifts in each seizure
    ax = findobj(h, 'type', 'axes');
    linkaxes(ax, 'y');
    ylim(ax, [0 max(dat_shifts.yy) + 1])
    box(ax, 'off')
    
    
end

function PP = parse_inputs_(varargin)
    
    PP = inputParser;
    PP.KeepUnmatched = true;
    
    p_ = @(varargin) PP.addParameter(varargin{:});
    p_('orientation', 'vertical', @(x) ismember(x, {'horizontal', 'vertical'}));
        % change figure orientation. <combine_stable_discrete_figs_>
    p_('show_merges', false);  % used in <merge_phases_>
    p_('create_fig', false);  % used in <plot_stable_modes_>
    p_('color_modes_by', 'phase_adj', @(x) ismember(x, {'phase', 'phase_adj', 'time_rel'}));
    p_('center_by', 'first', @(x) ismember(x, {'longest', 'iw', 'first', 'none'}));
    p_('cmap', 'metric', @(s) ismember(s, {'metric', 'red'}));
    p_('modes2show', 'all', @(s) ismember(s, {'all', 'onsets'}));
    p_('thresh_merge_phases', pi/4, @(x) x >= 0 && x <= pi); % 2° <adjust_phases_>
    p_('thresh_small_shifts', 0, @(x) x >= 0 && x <= pi); % 2° <summarize_data_>
    p_('fig_suffix', '');

    parse(PP, varargin{:});
    
    % Create a suffix for the figure to show the params used
    non_defaults = setdiff(PP.Parameters, PP.UsingDefaults);
    if ismember('fig_suffix', non_defaults)
        fig_suffix = PP.Results.fig_suffix;
    elseif isempty(non_defaults), fig_suffix = ''; 
    else
        fig_suffix = '';
        for ff = string(non_defaults(:))'
            val = PP.Results.(ff);
            if isnumeric(val) && isscalar(val)
                fig_suffix = sprintf('%s_%s_%0.3g', ...
                    fig_suffix, ff, PP.Results.(ff));
            elseif ischar(val) || isstring(val)
                fig_suffix = sprintf('%s_%s_%s', ...
                    fig_suffix, ff, PP.Results.(ff));
            elseif strcmpi(ff, 'create_fig') % do nothing
            else
                error('Undefined suffix behavior')
            end
        end
        
    end

    PP = PP.Results;
    PP.fig_suffix = strrep(fig_suffix, '.', 'p');

end
    