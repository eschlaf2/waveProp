classdef MetricComparison < handle
   methods
       function MC = MetricComparison(WP, metrics)
           if nargin < 1
               sz = BVNY.load_seizures;
               N = height(sz);
               fname = @(ii) sprintf('%s_Seizure%d_fits.mat', ...
                   sz.patient{ii}, sz.seizure(ii));
               MC = arrayfun(@(ii) MetricComparison(fname(ii)), 1:N);
               return
           end
           if nargin < 2, metrics = []; end
           if ischar(WP), WP = WaveProp.load('files', dir(WP)); end
           MC.Fit = WP;
           if ~isempty(metrics)
               MC.Metrics = metrics;
           end
       end
       
       function dz = get.diff(MC)
           dz = MC.Data1.diff(MC.Data2);
       end
       function mn = get.mean(MC)
           mn = circ_mean(MC.diff, [], [], 'omitnan');
       end
       function md = get.median(MC)
           alpha = MC.diff(~isnan(MC.diff));
           md = circ_median(alpha);
           md = angle(exp(1j * md));
       end
       function c = get.conf(MC)
           [mn, ul] = circ_mean(MC.diff, [], [], 'omitnan');
           c = ul - mn;
       end
       function c = get.conf_onesided(MC)
           theta = MC.diff;
           theta(isnan(theta)) = [];
           mn = abs(MC.mean);
           t = circ_confmean(theta,0.1);
           c = mn + t;
       end
       function sd = get.std(MC)
           sd = circ_std(MC.diff, [], [], 'omitnan');
       end
       function out = get.xcorr(MC)
           % xcorr = [xc, pval]
           data = MC.direction;
           [xc, pval] = circ_corrcc(data(:, 1), data(:, 2));
           out = [xc, pval];
       end
       function out = get.kuiper_test(MC)
            % kuiper_test = [pval, k, K]
            % [pval, k, K] = circ_kuipertest(data(:, 1), data(:, 2), 100, 1);
            [pval, k, K] = circ_kuipertest(MC.Data1.Direction, MC.Data2.Direction);
            out = [pval, k, K];
       end
       function dir = get.direction(MC)
           % returns non-nan directions
           data1 = MC.Data1;
           data2 = MC.Data2;
           if numel(data1.Direction) > numel(data2.Direction)
               data1 = data1.resample_t0(data2.time);
           else
               data2 = data2.resample_t0(data1.time);
           end
           
           dir = [data1.Direction, data2.Direction];
           
           
           dir(any(isnan(dir), 2), :) = [];
       end
       function out = get.mtest(MC)
           % mtest = [h mu ul ll]
           alpha = MC.diff;
           alpha(isnan(alpha)) = [];
           [h mu ul ll] = circ_mtest(alpha, 0);
           out = [h mu ul ll];
       end
       function N = get.proportion_valid(MC)
           N = numfinite(MC.diff) ./ numel(MC.diff);
       end
       function N = get.num_valid(MC)
           N = numfinite(MC.diff);
       end
       function pval = get.wwtest(MC)
           data1 = MC.Data1.Direction;
           data2 = MC.Data2.Direction;
           data1(isnan(data1)) = [];
           data2(isnan(data2)) = [];
           pval = circ_wwtest(data1, data2);
       end
       function kappa = get.kappa(MC)
           % kappa = circ_kappa(MC.diff)
           alpha = MC.diff;
           alpha(isnan(alpha)) = [];
           kappa = circ_kappa(alpha);
       end
       function pval = get.rtest(MC)
           alpha = MC.diff;
           alpha(isnan(alpha)) = [];
           pval = circ_rtest(alpha);
       end
       function out = get.ktest(MC)
           % ktest = [p, f, rbar]
           data = MC.direction;
           n = size(data, 1);
           R1 = n*circ_r(data(:, 1));
           R2 = n*circ_r(data(:, 2));
           rbar = (R1+R2)/(2*n);
           if rbar < 0.7, out = [nan nan rbar]; return, end
           [p, f] = circ_ktest(data(:, 1), data(:, 2));
           out = [p, f, rbar];
       end
       function data1 = get.Data1(MC)
           data1 = MC.Fit.(MC.M1);
       end
       function data2 = get.Data2(MC)
           data2 = MC.Fit.(MC.M2);
       end
       function p = get.Patient(MC)
           info = strsplit(MC.Fit.Name{:}, '_');
           p = info{1};
       end
       function s = get.Seizure(MC)
           info = strsplit(MC.Fit.Name{:}, '_');
           s = info(contains(lower(info), 'seizure'));
           s = str2double(s{:}(8:end));
       end
       function m1 = M1(MC), m1 = MC.Metrics{1}; end
       function m2 = M2(MC), m2 = MC.Metrics{2}; end
       function set_metrics(MC, value), for m = MC, m.Metrics = value; end, end
       function set(MC, field, value)
           ss = validatestring(field, properties(MC));
           for m = MC, m.(ss) = value; end
       end
       function agree = get.StrongAgree(MC)
           agree = abs(MC.mean) + MC.conf < MC.H0;
       end
       function agree = get.WeakAgree(MC)
           agree = ~MC.StrongAgree & ...
               abs(MC.mean) < MC.H0 & ...
               abs(MC.mean) + MC.conf < pi/2;
       end
       function agree = get.Agree(MC)
           agree = MC.StrongAgree | MC.WeakAgree;
       end
       function biased = get.Biased(MC)
           biased = ~MC.Agree & MC.conf < pi/3;
       end
       function disagree = get.Disagree(MC)
           disagree = ~MC.Agree & ~MC.Biased;
       end
       
        function xx = xlocs(F)
            patient = {F.Patient};
            G = findgroups(patient);
            dG = G(2:end) - G(1:end - 1);
            xx = (1:numel(patient)) + cumsum([0 logical(dG)]);  % add a gap between patients
        end
       
    function ax = plot(MC, varargin)
        % Plot the circular mean and standard deviation comparing chosen metrics
        % for wave direction. Highlight individual patients and seizures.
        % USAGE: 
        %    ax = plot(MC, ax=axes(figure), ::'noreorder'::, varargin)

        [MC, ax, line_color, args] = parse_inputs_(MC, varargin);

        % Extract variables from table
        theta = cat(1, MC.mean);
        cconf = cat(1, MC.conf);
        seizure = cat(1, MC.Seizure);

        nanconf = isnan(cconf);
        cconf(nanconf) = pi;

        % Compute lower and upper bounds (circular confidence)
        lowCI = theta - cconf;
        hiCI = theta + cconf;
        

        %% Make the figure

        % Plot the means first so it's easy to label
%         G = findgroups(patient);
%         X = (1:N) + cumsum([0 logical(diff(G))]);  % add a gap between patients
        X = MC.xlocs;
        yy = nan(1, max(X));  % gaps should have nan values
        yy(X) = theta;

        plot(ax, yy, 'color', line_color, args{:}, 'tag', 'main');  % plot connecting lines
        set(ax,'TickLength',[0 0], ...
            'XTick', X, ...
            'XTickLabel', seizure, ...
            'YGrid','on', ...  
            'YTick', (-pi/4:pi/4:pi/4), ...
            'YTickLabel', [], ...
            'box', 'off', ...
            'nextplot', 'add', ...
            'xlim', [0 max(X)+1], ...
            'ylim', [-pi pi]);

        include_zero = lowCI <= 0 & hiCI >=0 & ~nanconf;
        theta_small = abs(theta) < pi/4 & ~nanconf;
%         disagree = abs(theta)+cconf > pi/4 | nanconf;
        disagree = cat(1, MC.Disagree) | cat(1, MC.Biased);

        x = [1;1] * X;
        y = [lowCI'; hiCI'];
        toohi = y(2, :) > pi;
        toolow = y(1, :) < -pi;

        plot(ax, x, y, 'color', line_color, 'tag', 'CIs');  % plot CIs
        plot(ax, x(:, toohi), ...  % Correct wrap-arounds
            [0*y(1, toohi) - pi; y(2, toohi) - 2*pi], ...
            'color', line_color, 'tag', 'CIs');
        plot(ax, x(:, toolow), ...
            [0*y(2, toolow) + pi; y(1, toolow) + 2*pi], ...
            'color', line_color, 'tag', 'CIs'); 
        scatter(ax, X, theta, 50, line_color, 'filled', 'tag', 'main');  % Put a dot at each mean
        scatter(ax, ...  % Show nans as gray diamonds
            X(nanconf), theta(nanconf), 90, .4*[1 1 1], ...
            'd', 'filled', 'tag', 'nanconf');
        scatter(ax, ...  % indicate where difference includes zero with a circle
            X(include_zero), theta(include_zero), 90, [1 0 0], 'tag', 'CI_incl_0');  
        plot(ax, ...    % indicate small differences with a *
            X(theta_small), theta(theta_small), 'r*', 'tag', 'small_diff');
        scatter(ax, ...    % indicate disagree differences with a star
            X(disagree), theta(disagree)*0 + .75*pi, 25, [0 0 0], ...
            'filled', 'Marker', 'p', 'tag', 'disagree');
        ylabel(ax, sprintf('%s - %s', ...
            MC(1).M2, MC(1).M1 ...
            ));
        hold(ax, 'off');
        yline(ax, 0);


%         annotation('textbox', ...
%             'string', {'\circ := CI includes zero'; '* := \theta < \pi/4'}, ...
%             'Position', [.01 .9 .05 .04], 'FitBoxToText', true);
        ax.Tag = [MC(1).M1 '_' MC(1).M2];

    end
    function color = get.Color(self)
        mtc = self.Metrics;
        cmap = lines;
        %         cmap = dark2; 
        testfun =@(x) all(contains(mtc, x));
        if testfun({'M', 'E'})
            color = cmap(1, :);
        elseif testfun({'M', 'D1xwh'})
            color = cmap(2, :);
        elseif testfun({'E', 'D1xwh'})
            color = cmap(3, :);
        else
            color = cmap(4, :);
        end
    end

       
   end
   properties
       Fit 
       Metrics (1, 2) = {'M', 'E'}
       Color
       H0 = pi/4
   end
   properties (Dependent = true, SetAccess = private)
       Data1
       Data2
       Patient
       Seizure
       Agree
       Disagree
       WeakAgree
       StrongAgree
       Biased

       diff
       std
       mean
       median
       conf
       xcorr
       ktest
       direction
       conf_onesided  % one-sided test (mean < theta, p < .05)
       kuiper_test  % analogue of KS test
       mtest  % H0: mean of differences == 0
       wwtest % H0: the two populations have equal means (returns pval)
       kappa  % estimate of concentration parameter of VM dist
       rtest  % Rayleigh test for uniformity. H0: the population is uniformly distributed about the circle
       proportion_valid  % proportion of non-nan (finite) comparisons
       num_valid  % number of non-nan (finite) comparisons
   end
end

%% LOCAL functions


function [MC, ax, line_color, args] = parse_inputs_(MC, args)
    line_color = lines(1);
    ax = [];
    isax = cellfun(@(a) isa(a, 'matlab.graphics.axis.Axes'), args);
    if any(isax), ax = args{isax}; args(isax) = []; end
    if isempty(ax), ax = axes(figure); end
    
    char_args = cellfun(@ischar, args);
    noreorderA = strcmpi(args(char_args), 'noreorder');
    if any(noreorderA)
        args(char_args(noreorderA)) = []; 
    else  % reorder MC
        patients = lower({MC.Patient});
        mask = ismember(patients, {'sim', 'scm', 'fs', 'sw'});
        MC = [MC(~mask), MC(mask)];
    end
    
    for ii = 1:2:length(args)  % find color NV pair
       if strcmpi(args{ii}, 'color'), line_color = args{ii + 1}; break; end 
    end
    
end


