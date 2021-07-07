function [ax, summ] = stats_summary(F, ff)

ksargs = F.KSArgs;
ax = axes(figure);
metrics = struct2cell(F.Metrics)';
metrics(cellfun(@isempty, metrics)) = [];   % no longer using 'E' or 'BOS*' methods (4/13/21)
fits = F.Fits;

use_circ_ks = {'Phi_mean_early', 'Phi_mean_late', 'Phi_early', 'Phi_late'};

summ = struct;
for mm = string(metrics)
    if contains(mm, 'dt'), ls = '--'; else, ls = '-'; end
    switch ff
        case "Phi_early"
            fit = WaveProp.resize_obj(fits.(mm));
            data = fit.Direction(fit.early);
        case "Phi_late"
            fit = WaveProp.resize_obj(fits.(mm));
            data = fit.Direction(fit.late);
        otherwise
            for ii = length(fits.Name):-1:1
                data(ii) = fits.(mm)(ii).(ff);
            end
    end
    summ.(mm) = data;
    if all(isnan(data)), continue; end
    
    
    try
        if isempty(data), fprintf('%s, %s bw: empty\n', ff, mm); continue; end
        if ismember(ff, use_circ_ks)
            [d, xi, bw] = circ_ksdens(data, ksargs.(ff){:});
            xi = xi/pi * 180;
            d = d/pi * 180;
        else
            [d, xi, bw] = ksdensity(data, ksargs.(ff){:});
        end
        
        fprintf('%s, %s bw: %.04f\n', ff, mm, bw);
        plot(ax, xi, d, F.Style.(mm));
        hold(ax, 'on');
    catch ME
        if ~strcmpi(ME.identifier, 'stats:mvksdensity:NanX')
            rethrow(ME)
        end
    end 
end     
            
hold(ax, 'off')
legend(ax)
title(ax, strrep(ff, '_', ' '))
ylabel('PDF')
switch ff
    case "N_detections_early"
        xlabel(ax, 'Count');
    case use_circ_ks
        xlabel(ax, 'Direction (°)');
        if strcmpi(ff, 'phi_mean_early')
            xticks(ax, (-180:90:180));
        end
    case "First_detection"
        xlabel(ax, 'Time (s)');
    case "Phi_std_late"
        ln = ax.Children;
        set(ln, 'xdata', ln(1).XData / pi * 180);
        xlabel(ax, 'St. Dev. (°)');
end

ax.Tag = ff;
