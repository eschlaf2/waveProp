%
% plot time-frequency of a signal over time using wavelet decomposition
% 

function [wavelet] = plot_wavelet_tf (wavelet, freqs, GRIDSIZE, yticks, cbar)

if nargin < 1
    warning ('Usage:  plot_wavelet_tf (wavelet, [freq vector], [gridsize], [yticks ([]=off)])');
    return;
end
if nargin < 2 || isempty(freqs)
    freqs = [min(wavelet.frequencies) max(wavelet.frequencies)];
end

if nargin < 3 || isempty(GRIDSIZE)
    GRIDSIZE = 400;
end

axis_off = 0;
if nargin < 4
    yticks = freqs;
else
    if isempty(yticks)
        axis_off = 1;
    end   
end
if nargin < 5 || isempty(cbar)
    cbar = 0;
end

% set Y resolution to nearest integer multiple of gridsize
nf = length(wavelet.frequencies);
yscale = max (1, floor(GRIDSIZE/nf));
YRESOLUTION = nf*yscale;

% find the y position corresponding to the ticks in the resampled data

display_freqs_i = find(wavelet.frequencies >= freqs(1) & wavelet.frequencies <= freqs(2));
display_freqs = wavelet.frequencies(display_freqs_i);

for t = 1:length(yticks)
    y = find(yticks(t) <= display_freqs, 1, 'first');
    if isempty(y)
        if t == 1, y = 1; end;
        if t > 1, y = length(yticks); end
    end
    display_ticks(t) = y*yscale;
    display_tick_labels{t} = num2str(yticks(t));
end    

amplitude = wavelet.corrected_amplitude(display_freqs_i,:);
nsamples = size(amplitude,2);

if nsamples < GRIDSIZE
    % apply cubic spline interpolation
    clear x y; [x,y] = meshgrid (1:size(amplitude,2),1:size(amplitude,1));
    clear xi yi; [xi, yi] = meshgrid (linspace(1,size(amplitude,2),GRIDSIZE), ...
        linspace(1,size(amplitude,1),YRESOLUTION));
    clear zi; zi = interp2 (x, y, amplitude, xi, yi, 'cubic');
else
    % downsample the X resolution but have to do it one channel at a time
    for k = 1:size(amplitude,1)
       resampled_amplitude(k,:) = resample(amplitude(k,:), GRIDSIZE, nsamples);
    end
    clear x y; [x,y] = meshgrid (1:size(resampled_amplitude,2),1:size(amplitude,1));
    clear xi yi; [xi, yi] = meshgrid (linspace(1,size(resampled_amplitude,2),GRIDSIZE), ...
        linspace(1,size(amplitude,1),YRESOLUTION));
    clear zi; zi = interp2 (x, y, resampled_amplitude, xi, yi, 'cubic');
end
    
nsecs = length(wavelet.wave)/wavelet.samp_rate;
gridlength = length(zi);
xticks = [1:gridlength]./wavelet.samp_rate;
if nsamples > GRIDSIZE
    xticks = xticks.*(nsamples/GRIDSIZE);
end

ax = axes(figure('color','w')); 
imagesc(ax, xticks, display_freqs, zi);
set(ax,'YDir','normal', 'YScale', 'log');
% hold on;
if 1%cbar > 0
    minval = min(nonzeros(zi));
    maxval = 0.05*max(nonzeros(zi));
    caxis(ax, [minval maxval]);
    colorbar(ax);
end
if axis_off
    axis('off');
else
%     set(ax,'YTick',yticks);
%     set(ax,'YTickLabel',display_tick_labels,'FontName','Helvetica','Fontsize',10);
	set(ax,'FontName','Helvetica','Fontsize',10);
%     set(ax,'xtick',[]);
end
% hold on;
% in order to do xticks, need to make up a second set of axes
% ax1 = gca;
% ax2 = axes('position',get(ax1,'position'),'color','none');
% hold on;
% plot(ax2,xticks,zeros(1,length(xticks)),'.');
% v = axis;
% set(gca,'ytick',[]);
% set(gca,'ylim',[-2 -1]);
% set(gca,'xlim',[v(1) v(2)]);



