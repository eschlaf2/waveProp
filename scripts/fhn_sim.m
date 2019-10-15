defvar(who, 'SAVE', false);
mx = zeros(50e4, 1);  % temporary - watch max values

ncols=150;                               % Number of columns in domain
nrows=150;                               % Number of rows in domain
dt=.5;                                   % Time step (ms)
dur=80e3/dt;                             % Number of time steps
mindur=45e3/dt;
% h=2.0;                                   % Grid size
% h2=h^2;
Iex=.5;                                  % Amplitude of external current
% mu=1.0;                                  % Anisotropy
G = .7;                                  % Conductance
% Gx=.7; Gy=G/mu;                          % Conductances
a = 0.45; b = 0.8; phi = 0.15;  % FHN model parameters
dvdt = @(v, r, I) v - v.^3/3 - r + I;  % Fast dynamics
drdt = @(v, r) phi * (v + a - b * r);  % Slow dynamics;

v=zeros(nrows,ncols);                    % Initialize voltage array
r=v;                                     % Initialize refractoriness array
skipfactor=50;                           % Save every nth frame
chunk = 5 * 1e3/dt;  % Save results in chunks

fname = strrep(sprintf('spiral_wave_%d', StimProtocol), '-', 'M');
basename = sprintf('%s%s%s', fname, filesep, fname);
res = nan(nrows, ncols, chunk, 'single');

if SAVE
	v_out=zeros(nrows,ncols, dur, 'single'); 
	r_out = v_out;
end
	
if ~exist(fname, 'dir'), mkdir(fname), end

% Set initial stim current and pattern
iex=zeros(nrows,ncols);

% Define the seizure (stim location and duration)
t = (dt / 1e3:dt / 1e3:80)';  % time (s)
seizure = define_seizure(t);
mea = floor([(nrows/2) (ncols/2)]);

% switch StimProtocol
% 	case {-1, 1, 3}
% % 		iex(62:67,62:67)=Iex;
% 		iex(1:6,1:6) = Iex;
% 	case {0, 2}
% 		iex(:,1)=Iex;
% end

% Setup image
ih=imagesc(v); set(ih,'cdatamapping','direct')
colormap(hot); axis image off; th=title('');
set(gcf,'position',[500 600 256 256],'color',[1 1 1],'menubar','none')
mov(dur/skipfactor) = getframe;

n=0;                 % Counter for time loop
n = 31e3/dt;
k=0;                 % Counter for movie frames
done=0;              % Flag for while loop


while ~done          % Time loop
    n = n + 1;
	inds = round(mea + seizure(n, [1 2]));
	iex(inds(1) + (0:5), inds(2) + (0:5)) = Iex * seizure(n, 3);
    
    % Create padded v matrix to incorporate Newman boundary conditions 
	vv = padarray(v(2:end-1, 2:end-1), [1 1], 'replicate', 'both');
    
    % Update v
	L = del2(vv);  % Laplacian
	v_new=v + dvdt(v, r, iex+G*L)*dt;
	
	mx(n) = max(v, [], 'all');
	if any(v_new(:) > 1e3)  % decrease dt when growth is too steep
% 		error('High v');
		[t, vn] = ode45(...
			@(t, y) dvdt(y, r(:), iex(:)+G*L(:)), ...
			[0 dt], ...  % (t range)
			v(:)); % y0
		v_new = reshape(vn(end, :), size(v_new));
	end
	
    % Update r
    r=r + drdt(v, r)*dt;
	v = v_new;
	
	% Store r and v
	if SAVE
		r_out(:, :, n)=r;
		v_out(:, :, n)=v; 
	end
    
    % Map voltage to image grayscale value
    m=1+round(63*v); m=max(m,1); m=min(m,64);
    
    % Update image and text 
    
    % Write every 500th frame to movie 
    if rem(n,skipfactor)==0
        set(ih,'cdata',m);
        set(th,'string',sprintf('%d  %0.2f   %0.2f',n,v(1,1),r(1,1)))
        drawnow
    
		res(:, :, mod(k, chunk) + 1) = v;
		if rem(k, chunk)==0, save(sprintf('%s_%06d', basename, k), 'res', 't'); end
		k = k + 1;
%         mov(k)=getframe;

% 		imwrite(uint16(v * 1e3), sprintf('%s_%06d.png', basename, k))
% 		print(gcf, sprintf('%s%s%s_%d', fname, filesep, fname, n), '-dpng');
		
    end
    
    done=(n > dur);
    if n > mindur && max(v(:)) < 1.0e-4, done=1; end      % If activation extinguishes, quit early.
    if ~isempty(get(gcf,'userdata')), done=1; end % Quit if user clicks on 'Quit' button.
end

res = res(:, :, 1:k);
if SAVE
	v_out = v_out(:, :, 1:n);
	r_out = r_out(:, :, 1:n);
	save(basename, 'res', 'v_out', 'r_out');
end

% % Write movie as AVI
% mov = mov(1:k-1);
% vid = VideoWriter([basename '.avi']);
% open(vid);
% writeVideo(vid, mov);
% close(vid);

% if isunix, sep='/'; else sep='\'; end
% [fn,pn]=uiputfile([pwd sep 'SpiralWaves.avi'],'Save movie as:');
% if ischar(fn)
%     movie2avi(mov,[pn fn],'quality',75)
% else
%     disp('User pressed cancel')
% end

close(gcf)

function seizure = define_seizure(t)

tau = [10 20 20 10 10]';  % Seizure stages
tau = cumsum(tau);

% Linear mapping function [t0, tf] -> [0, 1]
T =@(t0, tf, t) (t - t0) ./ (tf - t0);  

% Stim duration function
dss = .02;  % Single spike
d0 = .07;  % Short burst
df = .2;  % Long burst
D =@(t) ...  
	dss * (t > tau(1) & t <= tau(3)) + ...
	d0 * (t > tau(3) & t <= tau(4)) + ...
	(d0 + (df - d0) * T(tau(4), tau(5), t)) .* (t > tau(4) & t <= tau(5));

% Radius function
r0 = 10;
rf = 70;
R =@(t) ...
	(r0 + (rf - r0) * T(tau(1), tau(2), t)) .* (t > tau(1) & t <= tau(2)) + ...
	(rf) * (t > tau(2));

% Angle function
theta0 = -pi/2;
thetaf = 0;
theta =@(t) ...
	theta0 * (t <= tau(2)) + ...
	(theta0 + (thetaf - theta0) .* T(tau(2), tau(3), t)) .* (t > tau(2) & t <= tau(3)) + ...
	thetaf * (t > tau(3));

% Stim times
N = 100;  % Number of discharges
noise = .01;  % sd (in ms)
Ss = linspace(tau(1), tau(5), N)' + noise*randn(N, 1);  % starts
Sf = Ss + D(Ss);  % ... and ends

% Prep output
posx = R(t) .* cos(theta(t));  % Stim position (relative to mea; pos in units of electrodes)
posy = R(t) .* sin(theta(t));
I = mod(sum((t >= Ss') + (t <= Sf'), 2), 2);  % Stim indicator

% Stimulus map
seizure = [posx posy I];
end