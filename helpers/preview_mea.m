function preview_mea(mea, t_start, skipfactor, pause_length)
% Inputs:
%    mea (required)
%    t_start (default: 0): time (s) to start preview
%    skipfactor (default: 100): how to skip through frames
%    pause_length (default: 0): pause time (s) between rendering frames

defvar('skipfactor', 100);
defvar('pause_length', 0);
defvar('t_start', 0);

figure();
im = nan(max(mea.Position));
inds = sub2ind(size(im), mea.Position(:, 1), mea.Position(:, 2));
im(inds) = mea.Data(1, :);
ih = imagesc(im);
colormap bone;
th = title('');
t = mea.Time();

nt = length(mea.Data);
for ii = find(t >= t_start, 1):skipfactor:nt
	im(inds) = mea.Data(ii, :);
	set(ih, 'cdata', im); set(th, 'string', num2str(t(ii)));
	drawnow;
	pause(pause_length)
end
