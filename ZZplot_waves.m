function [] = plot_wave_fit(position, data, beta)

X = position(:, 1);
Y = position(:, 2);

p1 = subplot(1,3,1:2);
scatter3(p1, X, Y, data, 200, data, 'filled');hold on;        

[XX, YY] = meshgrid(sort(unique(X)),sort(unique(Y)));
Z = double(position * beta(1:2) + beta(end));
f = scatteredInterpolant(X, Y, Z);
ZZ = f(XX,YY);
mesh(p1, XX, YY, ZZ, ...
	'FaceColor', 'interp', 'FaceAlpha', 0.8, 'LineStyle', 'none') %interpolated

legend(p1, 'Data','Regression fit');
xlabel('cm');ylabel('cm');zlabel('Second');
hold off;
% Plot the projection along the velocity axis
p2 = subplot(133);
P_v_axis = position*beta(1:2)/norm(beta(1:2));
plot(p2, P_v_axis, data, '.');
hold on;
plot(p2, P_v_axis, Z);
title('Projection along the velocity vector');
xlabel('cm');
ylabel('Second'); colormap(p2, 'hot')
hold off;

% % Plot shuffled residual error distribution if L1-regression is used
% if strcmpi(p.Results.Lossfun,'L1')
% 	figure;
% 	hist(E_s);hold on;
% 	plot([E E],get(gca,'YLim'),'r');hold off
% 	title(['Mean residual error, p = ' num2str(pValue)]);
% end  

end 
