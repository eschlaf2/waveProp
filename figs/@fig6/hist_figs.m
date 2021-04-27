function hist_figs(F, pat, style)
    
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
    
    
    
%     iw_info = get_iw_info_(F.get_file(ff).name); 
    iw_info = get_iw_info_(F.get_file(ff));
%     if isempty(iw_info) || all(isnan(iw_info)), continue; end
    % center directions around IW
    if ~isnan(iw_info(end, 4))
        rotate_by = iw_info(end, 4);  % end to center MG49 around the second wave
        iw_info(:, 4) = angle(exp(1j.*(iw_info(:, 4) - rotate_by)));
    end
    
    ax_summ = ax(r, end);
    ax_summ.Tag = 'hist1d';  % used to be 'ksdens'
    set(ax(r, 1:end-1), 'tag', 'hist2D', 'nextplot', 'replacechildren');
    
    for c = 1:C-1
        mm = F.Metrics{c};
        W.(mm).RotateBy = rotate_by;
        W.(mm).MinFinite = BVNY.MinFinite;
        
        % make the original hist plots
        fit = W.(mm);
        fit.hist(F.Smoothing, false, false, h_temp); 
        
        % get info from the original main plot
        ax0 = findobj(h_temp, 'tag', 'Main');
        plt = ax0.Children;
        time = plt.XData;
        dir = plt.YData;
        tE = duration(ff);
%         mask = time >=0 & time <= tE;
        
        % update tmax
        if tE > tmax, tmax = tE; end
        
    
        % create the new 2D plot
        data = discretize_data_(plt);
        imagesc(ax(r, c), time, dir, data);
%         if COMBINE
%             rgb = 1-F.Style.(mm).color;
%             im0 = ind2rgb(data, map .* rgb);
% %             im0 = im0(:, mask, :);
%             if c == 1, im = im0; else, im = im + im0; end
%         else  % 'spread'
%             if r == 1
%                 title(ax(r, c), F.MetricNames.(mm), ...
%                     'color', F.Style.(mm).color); 
%             end
% %             spread_(plt);
% %             copyobj(plt, ax(r, c));
%             im = ind2rgb(data, 1-(map .* [1 1 1]));
%             imS = smooth_im_(im, A_RES);
%             imagesc(ax(r, c), time, dir, imS);
%         end
        
        show_iw_(ax(r, c), iw_info, F.Style.(mm).color);
        
        % make summary subplot
        update_summary_plot_(h_temp, ax_summ, F.Style.(mm))
    end
    

    
    
    

    if COMBINE
%         imagesc(ax(r, 1), time(mask), dir, transform_(im));
        imT = transform_(im);  % ** local fun **
        imS = smooth_im_(imT, A_RES);
        imagesc(ax(r, 1), time, dir, imS);  
%         imagesc(ax2, time, dir, imT);
    end
    set(ax(r, 1:end-1), 'xtick', [0 floor(tE)]);  
    

    
end     

prettify_(h, ax, tmax, P(ff))
h.Tag = [char(pat) '_' style];
close(h_temp)

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


function iw_info = get_iw_info_(fname)
    % iw_info = get_iw_info_(fname)
    % Output: [center ll ul phi]
    % Only show IWs in the first half of the seizure except in CUCX2_2
    fname0 = strsplit(fname, '_fits');
    load(['iw_mats/' fname0{1}], 'iw');
%     load('iw_info.mat', fname0);
    iw.reset;
    
    iw_info = nan(iw.num_waves, 5);
    for ii = 1:iw.num_waves
        
        if ii == iw.main_wave(iw.name)
            temp = BVNY.get_iw_info(iw.name);
            if isempty(temp), return; end
%             if temp(6) > .05, temp(4) = nan; end  % set phi to nan if fit is not significant
            iw_info(ii, :) = temp([1:4 6]);
            
        else
            iw.wave = ii;
            if ~iw.show, continue; end
            iw_info(ii, 1) = iw.center;
            iw_info(ii, 2:3) = iw.range;
            fit = iw.wave_fit;

            if fit.p < .05, phi = fit.phi; else, phi = nan; end
%             phi = iw.phi;
    %         if isnan(phi), phi = iw.wave_fit_alt.phi; end
            iw_info(ii, 4) = phi;
    %         iw_info(ii, 5) = iw.show;
        end
    end
    iw_info(isnan(iw_info(:, 1)), :) = [];
    if isempty(iw_info), iw_info = nan(iw.num_waves, 5); end
    
end


function update_summary_plot_(h_src, ax_dest, line_style)
    ax0 = findobj(h_src, 'tag', 'Summ');
    ln = ax0.Children;
%     ln.XData = rescale(ln.XData, 0, 1);
    set(ln, line_style);
    ln = copyobj(ln, ax_dest);
    grid(ax_dest, 'on');
    xticks(ax_dest, [0 1]);
    
    % Move D10 line to the back
    if strcmpi(ln.DisplayName, 'd10')
        ln.ZData = 0*ln.XData - 1;
    end
end


function prettify_(h, ax, tmax, patient) 
    COMBINE = size(ax, 2) == 2; 
    lbl = findobj(h, 'tag', 'label'); 
    lbl.YLabel.String = sprintf('Patient %d', patient); 
    
    set(findobj(h, 'tag', 'hist2D'), ... 
        'xlim', [-10 tmax+10], ... 
        'ticklength', [0 0], ... 
        'layer', 'top', ... 
        'xgrid', 'on', 'linewidth', 1);  
    set(ax, 'ytick', (-pi:pi/2:pi), 'yticklabel', [], ...
        'ygrid', 'on', 'ylim', 1.1*[-pi pi]);
    set(ax(:, 1), 'yticklabel', {'-\pi', '', '0', '', '\pi'}); 
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
%     time = plt.XData;
%     dataS = smoothdata(data, 2, 'gaussian', 3, 'SamplePoints', time);
    
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


