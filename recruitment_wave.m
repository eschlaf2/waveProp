function [frSm, ch, mea] = recruitment_wave(mea)
% Look at whether units get recruited to the seizure by a slowly moving
% wavefront.

[frSm, time, ch, mea] = smooth_firingRate(mea);                                  % Smooth the firing rate and limit to active channels
recruitmentInd = arrayfun(@(ii) find(frSm(:, ii) > 2, 1), 1:numel(ch));    % Find time points where the firing rate is 2sd above base
[recIndsSorted, so] = sort(recruitmentInd);                                % Get the order in which channels were recruited to seizure
position = get_position(mea, ch(so));                                      % Transform the position into coordinates


%% Plotting
figure(); fullwidth(true);
subplot(221);                                                              % Plot a 2D image of the recruitment times
data = time(recIndsSorted);
plot_data_v_location(position, data);

subplot(222);                                                              % Fit a wave to the recruitment times
[beta, ~, ~, ~, ~, p] = estimate_wave(data, position);
beta = circshift(beta, -1);
plot_wave_fit(position, data, beta, p);

subplot(2,2,3:4);                                                          % Plot the firing rate traces colored by order of recruitment
colorBy = recIndsSorted;
data = frSm(:, so);
plot_firing_rates(time, data, colorBy)
title([strrep(mea.Name, '_', ' ') ' Recruitment Wave'])


end

%% Local functions

function [] = plot_wave_fit(position, data, beta, p)
	if numel(beta) < 2, return, end
	X = position(:, 1);
	Y = position(:, 2);
	beta = beta(:);
	
	scatter3(X, Y, data, 200, data, 'filled');hold on;        

	[XX, YY] = meshgrid(sort(unique(X)),sort(unique(Y)));
	Z = double(position * beta(1:2) + beta(end));
	f = scatteredInterpolant(X, Y, Z);
	ZZ = f(XX,YY);
	mesh(XX, YY, ZZ, ...
		'FaceColor', 'interp', 'FaceAlpha', 0.8, 'LineStyle', 'none') %interpolated

	legend('Data','Regression fit');
	xlabel('X');ylabel('Y');zlabel('Second');
    title(sprintf('p = %0.2g', p))
	
	hold off;
end

function [] = plot_firing_rates(time, data, colorBy)
% Create a time series plot where traces are colored according to colorBy
% value
cmap = get(gcf, 'colormap');
cInds = round(interp1([min(colorBy), max(colorBy)], ...
    [1 length(cmap)], colorBy));
set(gca, 'colororder', cmap(cInds, :), 'nextplot', 'replacechildren')

if length(time) > length(data)
    time = interp1([1 length(data)], [time(1) time(end)], 1:length(data));
end
plot(time, data);
axis tight;
xlabel('Time [s]')
ylabel('Normalized firing rate');
end

function [] = plot_data_v_location(position, data)

img = nan(max(position));  % Initialize map to visualize recruitment timing
addy = sub2ind(size(img), position(:, 1), position(:, 2));  % assign a single index to each xy-location
img(addy) = data;  % Show when recruitment occured (in seconds)
imagesc(img')
axis xy; axis square;
xlabel('X');
ylabel('Y');
colorbar; 

end

function [position] = get_position(mea, ch)
% Get the position of coordinates of interest and transform it into
% coordinates (in case it is in some other unit)

P = mea.Position;
P(mea.BadChannels, :) = [];
Ps = P(ch, :);  % Limit to only coordinates of interest
if size(Ps, 1) == 1, position = (Ps - min(Ps)) + 1; return; end
position = (Ps - min(Ps)) ./ min(diff(unique(Ps))) + 1;  % Reassign coordinates
end

function [frSm, time, ch, mea] = smooth_firingRate(mea)
% Smooth the firing rate to highlight broad, slow changes in firing.
%     Downsample to 100 Hz
%     Limit to channels with mean firing rate at least 6 spikes per sec
%     Normalize

time = mea.Time; time = time();
try firingRate = mea.firingRate; 
catch; [firingRate, mea] = mua_firing_rate(mea);
end
ch = find(mean(firingRate) >= 6);  % Only use channels with mean firing rate at least 6 spikes per second
skipfactor = round(mea.SamplingRate / 100);
firingRate = downsample(firingRate, skipfactor);               % Downsample to 100 Hz
time = downsample(time, skipfactor);
frSm = smoothdata(firingRate(:, ch), 1, 'movmean', 1e3);                          % Smooth over a 10 second window;

mn = mean(frSm(time < 0, :));  % Compute the baseline mean
sd = std(frSm);  % ... and standard deviation from entire recording

frSm = (frSm - mn) ./ sd;  % Normalize to baseline
end
