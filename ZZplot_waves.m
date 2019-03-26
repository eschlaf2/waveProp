function [] = plot_wave_fit(position, data)

scatter3(position(:,1), position(:,2), 1000 * data, 200, 1000 * data, 'filled');
hold on    
plot3(position(center,1), position(center,2), 0, 'X');
x1fit = linspace(min(position(:,1)), max(position(:,1)), 10);
x2fit = linspace(min(position(:,2)), max(position(:,2)), 10);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = 1000*(b(1) + b(2)*X1FIT + b(3)*X2FIT);
mesh(X1FIT,X2FIT,YFIT, 'FaceColor', 'interp', 'FaceAlpha', 0.8, 'LineStyle', 'none')
xlabel('Electrode position (mm)');
ylabel('Electrodes position (mm)');
zlabel('Delay (ms)')
axis tight