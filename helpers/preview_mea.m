function preview_mea(mea, skipfactor, pause_length)

defvar('skipfactor', 100);
defvar('pause_length', 0);

figure();
im = nan(max(mea.Position));
inds = sub2ind(size(im), mea.Position(:, 1), mea.Position(:, 2));
im(inds) = mea.Data(1, :);
ih = imagesc(im);
th = title('');
t = mea.Time();

nt = length(mea.Data);
for ii = 1:skipfactor:nt
	im(inds) = mea.Data(ii, :);
	set(ih, 'cdata', im); set(th, 'string', num2str(t(ii)));
	drawnow;
	pause(pause_length)
end
