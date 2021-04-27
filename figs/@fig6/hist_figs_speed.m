function hist_figs_speed(F, pat, style)
    
pat = string(pat);
if nargin < 3 || isempty(style), style = 'spread'; end
style = validatestring(style, {'combined', 'spread'});
COMBINE = strcmpi(style, 'combined');
A_RES = 512;
which_file = find(strcmpi(F.Seizures.patient, pat));
R = numel(which_file);
C = numel(F.Metrics) + 1;
P = F.Seizures.patientAlt;
duration = F.Seizures.duration;

map = [0; linspace(0, 1, 128)'];


[h, ax] = F.create_fig(pat);
if COMBINE, ax = change_to_combined_fig_(ax); end  % ** local fun **

set(ax, 'colormap', 1-gray);
F.H.(pat) = h;
h_temp = figure;

tmax = 0;  % 140
for r = 1:R
    ff = which_file(r);
    W = WaveProp.load('files', F.get_file(ff), 'metrics', F.Metrics);
     
    iw_info = get_iw_info_(F.get_file(ff).name);
    
    ax_summ = ax(r, end);
    ax_summ.Tag = 'hist1d';  % used to be 'ksdens'
    set(ax(r, 1:end-1), 'tag', 'hist2D', 'nextplot', 'replacechildren');
    
    for c = 1:C-1
        mm = F.Metrics{c};
        
        % make the original hist plots
        fit = W.(mm);
        iw_info(:, 4) = nanmean(fit.Magnitude);
        hist_speed_(fit, F.Smoothing, h_temp);
        
        % get info from the original main plot
        ax0 = findobj(h_temp, 'tag', 'Main');
        plt = ax0.Children;
        time = plt.XData;
        speed = plt.YData;
        tE = duration(ff);
%         mask = time >=0 & time <= tE;
        
        % update tmax
        if tE > tmax, tmax = tE; end
        
    
        % create the new 2D plot
        if COMBINE
            rgb = 1-F.Style.(mm).color;
            data = discretize_data_(plt);  % ** local fun **
            im0 = ind2rgb(data, map .* rgb);
%             im0 = im0(:, mask, :);
            if c == 1, im = im0; else, im = im + im0; end
        else  % 'spread'
            if r == 1
                title(ax(r, c), F.MetricNames.(mm), ...
                    'color', F.Style.(mm).color); 
            end
%             spread_(plt);
%             copyobj(plt, ax(r, c));
            data = discretize_data_(plt);  % ** local fun **            
            im = ind2rgb(data, 1-(map .* [1 1 1]));
            imS = smooth_im_(im, A_RES);
            imagesc(ax(r, c), time, speed, imS);
        end
        
        show_iw_(ax(r, c), iw_info, F.Style.(mm).color);
        
        % make summary subplot
%         update_summary_plot_(h_temp, ax_summ, F.Style.(mm))
    end
    

    
    
    

    if COMBINE
%         imagesc(ax(r, 1), time(mask), dir, transform_(im));
        imT = transform_(im);  % ** local fun **
        imS = smooth_im_(imT, A_RES);
        imagesc(ax(r, 1), time, speed, imS);  
%         imagesc(ax2, time, dir, imT);
    end
    set(ax(r, 1:end-1), 'xtick', floor(tE));  
    

    
end     

prettify_(h, ax, tmax, P(ff))
h.Tag = [char(pat) '_' style];
close(h_temp)

end

%% local functions
function show_iw_(ax0, iw_info, color)
    
    assignin('base', 'iw_info', iw_info);
    N = size(iw_info, 1);
    sts = ax0.NextPlot;
    ax0.NextPlot = 'add';
    ylims = ax0.YLim;
    
    for ii = 1:N
        ln = xline(ax0, iw_info(ii, 1));  % IW center
        
        face_col = .5 * [1 1 1];
        [v1, v2] = ndgrid(iw_info(ii, 2:3), ylims);
        pp = patch(ax0, 'vertices', [v1(:) v2(:)], 'faces', [1 2 4 3], ...
            'facecolor', face_col, ...
            'facealpha', .3, ...
            'linestyle', 'none');
        dir = plot(ax0, iw_info(ii, 1), iw_info(ii, 4), '.', ...
            'color', color, ...
            'markersize', 15);
    end
    
    ax0.NextPlot = sts;
end


function iw_info = get_iw_info_(fname)
    % iw_info = get_iw_info_(fname)
    % Output: [center ll ul phi]
    % Only show IWs in the first half of the seizure except in CUCX2_2
    fname0 = strsplit(fname, '_fits');
    load(['iw_mats/' fname0{1}], 'iw');
    iw.reset;
    
    iw_info = nan(iw.num_waves, 4);
    for ii = 1:iw.num_waves
        iw.wave = ii;
        if iw.onset_rel > .5 && ~strcmpi(iw.name, 'CUCX2_Seizure2'), continue; end
        iw_info(ii, 1) = iw.center;
        iw_info(ii, 2:3) = iw.range;
        phi = iw.phi;
%         if isnan(phi), phi = iw.wave_fit_alt.phi; end
        iw_info(ii, 4) = phi;
%         iw_info(ii, 5) = iw.show;
    end
    iw_info(isnan(iw_info(:, 1)), :) = [];
    
end


function update_summary_plot_(h_src, ax_dest, line_style)
    ax0 = findobj(h_src, 'tag', 'Summ');
    ln = ax0.Children;
    ln.XData = rescale(ln.XData, 0, 1);
    set(ln, line_style);
    copyobj(ln, ax_dest);
    grid(ax_dest, 'on');
    xticks(ax_dest, [0 1]);
end


function prettify_(h, ax, tmax, patient) 
    COMBINE = size(ax, 2) == 2; 
    lbl = findobj(h, 'tag', 'label'); 
    lbl.YLabel.String = sprintf('Patient %d', patient); 
    linkaxes(ax, 'y')
    set(findobj(h, 'tag', 'hist2D'), ... 
        'xlim', [-10 tmax+10], ... 
        'ticklength', [0 0], ... 
        'layer', 'top', ... 
        'xgrid', 'on', 'linewidth', 1);  
    set(ax, 'ygrid', 'on');
    ymax = min(ax(1, 1).YLim(2), 2e3);
    set(ax, 'ytick', round(linspace(0, ymax/100, 3)) * 100, ...
        'yscale', 'log')
    set(ax(:, 2:end), 'yticklabel', '');
%     set(ax(:, 1), 'yticklabel', {'-\pi', '', '0', '', '\pi'}); 
    if COMBINE 
        lgd = legend(ax(1, end), 'location', 'east'); 
        lgd.Position(1) = .8; 
    end 
end 


function ax = change_to_combined_fig_(ax)
    for r = 1:size(ax, 1)
        ax(r, 1).Position(3) = ...
            sum(ax(r, 2).Position([1 3])) - ax(r, 1).Position(1);
        ax(r, end).Position(1) = ax(r, end - 1).Position(1);
    end
    delete(ax(:, 2:end-1));
    ax = ax(:, [1 end]);
end


function imS = smooth_im_(imT, A_RES)
    
    [Na, Nt, ~] = size(imT);
    [xx, yy, zz] = ndgrid(linspace(1, Na, A_RES), 1:Nt, 1:3);
    imS = interpn(imT, xx, yy, zz, 'spline');

end


function dataD = discretize_data_(plt)
    data = plt.ZData;
    
%     LL = linspace(CLIM(1), CLIM(2), 63);
    LL = quantile(data(:), linspace(.9, .999, 120));
    dataD = discretize(data, [-Inf, LL, Inf]);
    
end


function spread_(plt)
            
    plt.ZData(isnan(plt.ZData)) = 0;
%     plt.ZData = zscore(plt.ZData);
%     plt.LevelList = .5:.1:3;
    plt.LevelList = quantile(plt.ZData(:), linspace(.9, .999, 120));
    

end


function imT = transform_(im) 
%     H = -1/6;
    H = 0;
    S = 1.1;
    V = .9;
    
    imT = rgb2hsv(1 - im);
    
    if numel(imT) == 3
        imT(1) = mod(imT(1) + H, 1);
        imT(2) = S * imT(2);
        imT(3) = V * imT(3);
    else
        imT(:, :, 1) = mod(imT(:, :, 1) + H, 1);
        imT(:, :, 2) = S * imT(:, :, 2);
        
        val = imT(:, :, 3);
        mask = val < 1;
        val(mask) = V * val(mask);
        imT(:, :, 3) = val;
    end
    imT(imT < 0) = 0; imT(imT > 1) = 1;
    imT = hsv2rgb(imT);
end


function hist_speed_(fit, smoothing, h_temp)
    data = fit.Magnitude;
    t = fit.time;
    qq = quantile(data, [.01 .99]);
    xi = linspace(qq(1), qq(2), 200);
    
    t_sub = t(1):.1:t(end);
    N = numel(t_sub);
    
    H = zeros(numel(xi), N);
    
    for ii = 1:N
        inds = (t > t_sub(ii) - smoothing/2) & (t < t_sub(ii) + smoothing/2);
        dat = data(inds);
        if all(isnan(dat)), continue; end
        H(:, ii) = ksdensity(dat, xi);
    end
    
    set(0, 'currentfigure', h_temp);
    contourf(t_sub, xi, H);
    set(gca, 'tag', 'Main');
    
end

