function hist_figs2(F, pat)
    
pat = string(pat);

which_file = find(strcmpi(F.Seizures.patient, pat));
R = numel(which_file);
C = numel(F.Metrics) + 1;
P = F.Seizures.patientAlt;
duration = F.Seizures.duration;


[h, ax] = F.create_fig(pat);

F.H.(pat) = h;


tmax = 0;  % 140
for r = 1:R
    ff = which_file(r);
    W = WaveProp.load(F.get_file(ff), F.Metrics);
  
    % rotate some patients so that stable direction is not on the boundary
    % Started rotating by the IW angle rather than for path visibility.
    if 0 && ismember(lower(pat), {'c7', 'mg49', 'cucx4', 'cucx5', ...
            'sime', 'simj', 'simk', 'simm', 'simn'})
        rotate_by = pi;
    elseif strcmpi(pat, 'c5')
        rotate_by = pi/2;
    else
        rotate_by = 0;
    end
    
    
    
    [iw_info, main_wave] = get_iw_info_(F.get_file(ff));

    % center directions around IW
    if ~isnan(iw_info(main_wave, 4))
        rotate_by = iw_info(main_wave, 4);  % all directions should be set to nan except the main wave
        iw_info(:, 4) = angle(exp(1j.*(iw_info(:, 4) - rotate_by)));
    end
    
    ax_summ = ax(r, end);
    ax_summ.Tag = 'hist1d';  % used to be 'ksdens'
    ax_summ.NextPlot = 'add';
    set(ax(r, 1:end-1), 'tag', 'hist2D', 'nextplot', 'replacechildren');
    
    for c = 1:C-1
        mm = F.Metrics{c};
        W.(mm).RotateBy = rotate_by;
        W.(mm).MinFinite = BVNY.MinFinite(W.(mm));
        
        % make the direction rasters
        fit = W.(mm);
%         dir = fit.discharge_directions;
        dir = fit.Direction;
        aa = fit.direction_raster(dir, ax(r, c));
        aa.Tag = 'direction_raster';
        
        tE = duration(ff);
        
        % update tmax
        if tE > tmax, tmax = tE; end

        
%         show_iw_(ax(r, c), iw_info, F.Style.(mm).color);
        
        % make summary subplot
        update_summary_plot_(dir, ax_summ, F.Style.(mm))
    end
    

    set(ax(r, 1:end-1), 'xtick', [0 floor(tE)]);  
    

    
end     

prettify_(h, ax, tmax, pat)
h.Tag = char(pat);


end

%% local functions
function show_iw_(ax0, iw_info, color)
    
    if isempty(iw_info), return; end
    
    assignin('base', 'iw_info', iw_info);
    N = size(iw_info, 1);
    sts = ax0.NextPlot;
    ax0.NextPlot = 'add';
    ylims = ax0.YLim;
    
    for ii = 1:N
        if any(isnan(iw_info(ii, 1:3))), continue; end
        if iw_info(ii, 5) > 0.05, col = .5*[1 1 1]; else, col = color; end
        ln = xline(ax0, iw_info(ii, 1));  % IW center
        
        face_col = .5 * [1 1 1];
        [v1, v2] = ndgrid(iw_info(ii, 2:3), ylims);
        pp = patch(ax0, 'vertices', [v1(:) v2(:)], 'faces', [1 2 4 3], ...
            'facecolor', face_col, ...
            'facealpha', .3, ...
            'linestyle', 'none');
        
        dir = plot(ax0, iw_info(ii, 1), iw_info(ii, 4), '.', ...
            'color', col, ...
            'markersize', 15);
    end
    
    ax0.NextPlot = sts;
end

function [iw_info, main_wave] = get_iw_info_(fname)
    % iw_info = get_iw_info_(fname)
    % Output: [center ll ul phi phi_pval]
    fname0 = strsplit(fname, '_fits');
    [iw_stats, main_wave] = BVNY.get_iw_info(fname0{1});
    if isempty(iw_stats), iw_info = nan(1, 5); return; end
    n_waves = numel(iw_stats);
    iw_info = nan(n_waves, 5);
    
    empties = [];
    for ii = 1:n_waves
        temp = iw_stats{ii};
        if isempty(temp), empties = [empties; ii]; continue; end %#ok<AGROW>
        iw_info(ii, 1) = temp.onset;
        iw_info(ii, 2:3) = temp.range;        
        iw_info(ii, 4) = temp.wave_fit.phi; % this was in degrees at one point, but should be fixed now...
        iw_info(ii, 5) = temp.wave_fit.p;
    end
    
    iw_info(empties, :) = [];
end


function update_summary_plot_(dir, ax_dest, line_style)
    

    [f, xi] = circ_ksdens(dir);
    ln = plot(ax_dest, f, xi/pi*180);
    set(ln, line_style);
    grid(ax_dest, 'on');
    xticks(ax_dest, [0 1]);
    
end


function prettify_(h, ax, tmax, patient) 
    
    lbl = findobj(h, 'tag', 'label'); 
    lbl.YLabel.String = patient;
%     lbl.YLabel.String = sprintf('Patient %d', patient); 
    
    set(findobj(h, 'tag', 'direction_raster'), ... 
        'xlim', [-10 tmax+10], ... 
        'ticklength', [0 0], ... 
        'layer', 'top', ... 
        'xgrid', 'on', 'linewidth', 1);
    linkaxes(findobj(h, 'tag', 'direction_raster'), 'x');
    set(ax, 'ytick', -180:90:180, 'yticklabel', [], ...
        'ygrid', 'on', 'ylim', 1.1*[-180 180]);
    set(ax(:, 1), 'yticklabelmode', 'auto');
    
%     set(ax, 'ytick', (-pi:pi/2:pi), 'yticklabel', [], ...
%         'ygrid', 'on', 'ylim', 1.1*[-pi pi]);
%     set(ax(:, 1), 'yticklabel', {'-\pi', '', '0', '', '\pi'}); 
    
end 







