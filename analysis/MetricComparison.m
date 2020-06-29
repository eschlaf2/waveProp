classdef MetricComparison < handle
   methods
       function MC = MetricComparison(WP, metrics)
           if nargin < 1, return, end
           if nargin < 2, metrics = []; end
           if ischar(WP), WP = WaveProp.load('files', dir(WP)); end
           MC.Fit = WP;
           if ~isempty(metrics)
               MC.Metric1 = string(metrics(1));
               MC.Metric2 = string(metrics(2));
           end
       end
       
       function dz = get.diff(MC)
           dz = MC.Data1.diff(MC.Data2);
       end
       function mn = get.mean(MC)
           mn = circ_mean(MC.diff, [], [], 'omitnan');
       end
       function c = get.conf(MC)
           [mn, ul] = circ_mean(MC.diff, [], [], 'omitnan');
           c = ul - mn;
       end
       function sd = get.std(MC)
           sd = circ_std(MC.diff, [], [], 'omitnan');
       end
       
       function data1 = get.Data1(MC)
           data1 = MC.Fit.(MC.Metric1);
       end
       function data2 = get.Data2(MC)
           data2 = MC.Fit.(MC.Metric2);
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
       function set(MC, field, value)
           ss = validatestring(field, properties(MC));
           for m = MC, m.(ss) = value; end
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
        patient = {MC.Patient};
        seizure = cat(1, MC.Seizure);
        N = numel(MC);

        nanconf = isnan(cconf);
        cconf(nanconf) = pi;

        % Compute lower and upper bounds (circular confidence)
        lowCI = theta - cconf;
        hiCI = theta + cconf;


        %% Make the figure

        % Plot the means first so it's easy to label
        G = findgroups(patient);
        X = (1:N) + cumsum([0 logical(diff(G))]);  % add a gap between patients
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
        ylabel(ax, sprintf('%s - %s', ...
            MC(1).Metric2, MC(1).Metric1 ...
            ));
        hold(ax, 'off');
        yline(ax, 0);


%         annotation('textbox', ...
%             'string', {'\circ := CI includes zero'; '* := \theta < \pi/4'}, ...
%             'Position', [.01 .9 .05 .04], 'FitBoxToText', true);
        ax.Tag = [MC(1).Metric1 '_' MC(1).Metric2];

    end


       
   end
   properties
       Fit 
       Metric1
       Metric2
   end
   properties (Dependent = true, SetAccess = private)
       Data1
       Data2
       diff
       std
       mean
       conf
       Patient
       Seizure
   end
end

%% LOCAL functions


function [MC, ax, line_color, args] = parse_inputs_(MC, args)
    line_color = lines(1);
    
    isax = cellfun(@(a) isa(a, 'matlab.graphics.axis.Axes'), args);
    ax = args{isax};
    args(isax) = [];
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


