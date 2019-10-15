rng default
defvar(who, 'SAVE', false);

ncols=160;                               % Number of columns in domain
nrows=160;                               % Number of rows in domain
dt=.5;                                   % Time step (ms)
dur=80e3/dt;                             % Number of time steps
mindur=80e3/dt;
Iex=.5;                                  % Amplitude of external current
G = 1;                                  % Conductance
a = 0.48; b = 0.8; phi = 0.17;  % FHN model parameters
dvdt = @(v, r, I) v - v.^3/3 - r + I;  % Fast dynamics
drdt = @(v, r) phi * (v + a - b * r);  % Slow dynamics;

v=zeros(nrows,ncols);                    % Initialize voltage array
r=v;                                     % Initialize refractoriness array
skipfactor=20/dt;                           % Save every nth frame
chunk = 5 * 1e3/dt;  % Save results in chunks

fname = strrep(sprintf('spiral_wave_%d', StimProtocol), '-', 'M');
basename = sprintf('%s%s%s', fname, filesep, fname);


if SAVE
	v_out=zeros(nrows,ncols, chunk, 'int16'); 
	r_out = v_out;
end
	
if ~exist(fname, 'dir'), mkdir(fname), end

% Set initial stim current and pattern
iex=zeros(nrows,ncols);

% Define the seizure (stim location and duration)
t = (dt / 1e3:dt / 1e3:80)';  % time (s)
seizure = define_seizure(t);
mea = floor([(nrows/2) (ncols/2)]);
noise =@() randn(nrows,ncols) * .01;

% Setup image
ih=imagesc(v); set(ih,'cdatamapping','direct')
colormap(hot); axis image off; th=title('');
set(gcf,'position',[500 600 256 256],'color',[1 1 1],'menubar','none')
hold on;
patch('faces', [(1:4) 1], ...
	'vertices', [mea - [5 5]; mea + [5 -5]; mea + [5 5]; mea + [-5 5]], ...
	'edgecolor', [1 0 0], ...
	'facecolor', 'none', ...
	'linewidth', 2);
hold off;
mov(dur/skipfactor) = getframe;

n=0;                 % Counter for time loop
% n = 15e3/dt;
k=0;                 % Counter for movie frames
done=0;              % Flag for while loop


while ~done          % Time loop
    n = n + 1;
	inds = round(mea + seizure(n, [1 2]));
	iex(inds(1) + (-2:2), inds(2) + (-2:2)) = Iex * seizure(n, 3);
    
    % Create padded v matrix to incorporate Neumann boundary conditions
% 	vv=[[0 v(2,:) 0];[v(:,2) v v(:,end-1)];[0 v(end-1,:) 0]];
	vv = padarray(v(3:end-2, 3:end-2), [2 2], 'replicate', 'both');
    
    % Update v
	L = del2(vv);  % Laplacian
	v_new=v + dvdt(v, r, iex+G*L + noise())*dt;
	
% 	mx(n) = max(v, [], 'all');
	if any(v_new(:) > 1e3)  % decrease dt when growth is too steep
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
		v_out(:, :, mod(n-1, chunk)+1) = v*1e4;
		r_out(:, :, mod(n-1, chunk)+1) = r*1e4;
		if rem(n, chunk)==0
			fname = sprintf('%s_%06d', basename, n);
			fprintf('Saving %s ...', fname);
			save(fname, 'v_out', 'r_out');
			fprintf(' Done.\n')
		end
	end
    
    % Update image and text 
    if rem(n,skipfactor)==0
		% Map voltage to image grayscale value
		m=1+round(63*v); m=max(m,1); m=min(m,64);
        set(ih,'cdata',m);
        set(th,'string', ...
			sprintf('%0.1f  %0.2f   %0.2f',n*dt,v(mea(1),mea(2)),r(mea(1),mea(2))))
        drawnow
   		k = k + 1;
    end
    
    done=(n > dur);
    if n > mindur && max(v(:)) < 1.0e-4, done=1; end      % If activation extinguishes, quit early.
    if ~isempty(get(gcf,'userdata')), done=1; end % Quit if user clicks on 'Quit' button.
end

if SAVE
	v_out = v_out(:, :, 1:mod(n, chunk));
	r_out = r_out(:, :, 1:mod(n, chunk));
	save(sprintf('%s_%d', basename, n), 'v_out', 'r_out', 't');
end

close(gcf)

function seizure = define_seizure(t)

tau = [10 20 20 10 10]';  % Seizure stages
tau = cumsum(tau);

% Linear mapping function [t0, tf] -> [0, 1]
T =@(t0, tf, t) (t - t0) ./ (tf - t0);  

% Stim duration function
dss = .005;  % Single spike
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
N = 110;  % Number of discharges
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