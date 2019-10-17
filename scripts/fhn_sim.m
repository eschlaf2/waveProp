rng default
defvar(who, 'SAVE', false);

params = init_params();
struct2var(params.model);

v=zeros(nrows,ncols);                    % Initialize voltage array
r=v;                                     % Initialize refractoriness array

fname = ['FHN' filesep 'FHN'];
if ~exist('sim_num', 'var')
	sim_num = numel(dir(sprintf('%s_Seizure*mat', fname)));
end
if ischar(sim_num), sim_num = str2double(sim_num); end
basename = sprintf('%s%sFHN_%d', fname, filesep, sim_num);

if SAVE
	v_out=zeros(nrows,ncols, chunk, 'int16'); 
	r_out = v_out;
end
	
if ~exist(fname, 'dir'), mkdir(fname), end

% Set initial stim current and pattern
iex=zeros(nrows,ncols);

% Define the seizure (stim location and duration)
t = (dt/1e3 : dt/1e3 : 80)' - 10;  % time (s)
seizure = define_seizure(t, params.seizure);
mea_ctr = floor([(nrows/2) (ncols/2)]);
noise =@() randn(nrows,ncols) * .01;

% Setup image
ih=imagesc(v); set(ih,'cdatamapping','direct')
colormap(hot); axis image off; th=title('T');
set(gcf,'position',[500 600 256 256],'color',[1 1 1],'menubar','none')
hold on;
patch('faces', [(1:4) 1], ...
	'vertices', [mea_ctr - [5 5]; mea_ctr + [4 -5]; mea_ctr + [4 4]; mea_ctr + [-5 4]], ...
	'edgecolor', [1 0 0], ...
	'facecolor', 'none', ...
	'linewidth', 2);
hold off;
mov(dur/skipfactor) = getframe;

n=1;                 % Counter for time loop
% n = 24e3/dt;
k=1;                 % Counter for movie frames
done=0;              % Flag for while loop


while ~done          % Time loop
	inds = round(mea_ctr + seizure(n, [1 2]));
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
			sprintf('%0.1f  %0.2f   %0.2f',n*dt,v(mea_ctr(1),mea_ctr(2)),r(mea_ctr(1),mea_ctr(2))))
        drawnow
		mov(k) = getframe(ih.Parent.Parent);
   		k = k + 1;
    end
    
    done=(n >= dur);
    if n > mindur && max(v(:)) < 1.0e-4, done=1; end      % If activation extinguishes, quit early.
    if ~isempty(get(gcf,'userdata')), done=1; end % Quit if user clicks on 'Quit' button.
	n = n + 1;
end
n = n-1;

if SAVE
	v_out = v_out(:, :, 1:mod(n-1, chunk)+1);
	r_out = r_out(:, :, 1:mod(n-1, chunk)+1);
	save(sprintf('%s_%06d', basename, n), 'v_out', 'r_out', 't', 'mov', 'params');
end

close(gcf)

if SAVE
	mea = convert_to_mea_data(basename);
	mea.seizure = seizure;
	mea.fhn_params = params;
	save(mea.Path, '-struct', 'mea');
	convert_to_wavefit_data(mea);
end

function [params, p] = init_params(varargin)
p = inputParser;
add =@(varargin) addParameter(p, varargin{:});
p.KeepUnmatched = true;

% Integration parameters
add('ncols', 160);                               % Number of columns in domain
add('nrows', 160);                               % Number of rows in domain
add('dt', .5);                                   % Time step (ms)
add('dur', 160e3);                             % Number of time steps
add('mindur', 160e3);

% Coupling parameters
add('Iex', .25);                                  % Amplitude of external current
add('G', 1);                                  % Conductance

% FHN parameters
add('a', 0.48); 
add('b', 0.8); 
add('phi', 0.17);

% FHN dynamics
add('dvdt', @(v, r, I) v - v.^3/3 - r + I);  % Fast dynamics
add('drdt', '');

% Video and output params
add('skipfactor', 40);    % Save every nth frame
add('chunk', 1e4);  % Save results in chunks

parse(p, varargin{:});
params.model = p.Results;

if isempty(params.model.drdt)
	params.model.drdt = @(v, r) params.model.phi * (v + params.model.a - params.model.b * r);
end


%% Seizure params
p = inputParser;
add =@(varargin) addParameter(p, varargin{:});
p.KeepUnmatched = true;

% Length of each stage (in seconds)
validate_tau = @(x) numel(x) == 5 || numel(x) == 4;
add('tau', [0 20 20 10 10]', validate_tau);  % Seizure stages

% Stim duration params
add('dss', .005);  % Single spike
add('d0', .07);  % Short burst
add('df', .2);  % Long burst

% Radius function
add('r0', 10);
add('rf', 70);

% Angle function
add('theta0', -pi/2);
add('thetaf', 0);

% Stim times
add('N', 100);  % Number of discharges
add('noise', .01);  % sd (in ms)

parse(p, varargin{:});
params.seizure = p.Results;

end

function [seizure, R, theta] = define_seizure(t, params)

stages = {'pre', 'moving away', 'moving around', 'constant burst', 'increasing burst'};

struct2var(params);

tau = cumsum(tau(:));

% Linear mapping function [t0, tf] -> [0, 1]
T =@(t0, tf, t) (t - t0) ./ (tf - t0);  

% Stim duration function
D =@(t) ...  
	dss * (t > tau(1) & t <= tau(3)) + ...
	d0 * (t > tau(3) & t <= tau(4)) + ...
	(d0 + (df - d0) * T(tau(4), tau(5), t)) .* (t > tau(4) & t <= tau(5));

% Radius function
R =@(t) ...
	(r0 + (rf - r0) * T(tau(1), tau(2), t)) .* (t > tau(1) & t <= tau(2)) + ...
	(rf) * (t > tau(2));

% Angle function
theta =@(t) ...
	theta0 * (t <= tau(2)) + ...
	(theta0 + (thetaf - theta0) .* T(tau(2), tau(3), t)) .* (t > tau(2) & t <= tau(3)) + ...
	thetaf * (t > tau(3));

% Stim times
Ss = linspace(tau(1), tau(5), N)' + noise*randn(N, 1);  % starts
Sf = Ss + D(Ss);  % ... and ends

% Prep output
posx = R(t) .* cos(theta(t));  % Stim position (relative to mea; pos in units of electrodes)
posy = R(t) .* sin(theta(t));
Istim = mod(sum((t >= Ss') + (t <= Sf'), 2), 2);  % Stim indicator

% Stimulus map
seizure = [posx posy Istim];

end

function mea = convert_to_mea_data(fname, sim_num, t)
mea.BadChannels = [];
mea.Time = t;
mea.Name = ['FHN ' num2str(sim_num)];
mea.Padding = [10 10];
mea.Duration = t(end) - mea.Padding(2);
mea.SamplingRate = 2e3;
[xx, yy] = meshgrid(1:10, 1:10);
mea.Position = [xx(:) yy(:)];
mea.Path = sprintf(...
	'%s%s%s_Seizure%02d_Neuroport_10_10.mat', ...
	pwd, filesep, fname, sim_num);

files = dir([fname filesep '*mat']);
addpath(fname)

nf = length(files);

mea_addy = [80, 80];

win = 1./(sqrt((-3:3).^2 + (-3:3)'.^2) + .5);
win = win ./ sum(win(:));

for ii = 0:nf - 1
		load(files(ii+1).name, 'v_out');
		chunksize = size(v_out, 3);
		v_out = v_out(mea_addy(1) + (-8:7), ...
			mea_addy(2) + (-8:7), :);
		for jj = chunksize:-1:1
			temp(:, :, jj) = conv2(v_out(:, :, jj), win, 'valid');
		end
		
		mea.Data(:, :, (ii*chunksize + 1):((ii+1)*chunksize)) = -temp;
end
mea.Data = reshape(mea.Data, 100, [])';

end

function convert_to_wavefit_data(mea)

[~, fname, ~] = fileparts(mea.Path);
finfo = strsplit(fname, '_');
pat = finfo{1}; seizure = finfo{2}(8:end);
fits = matfile([fname '_wave_prop_1.mat'], 'writable', true);
plotnum = 0; wave_dir_polarplots;  % this should produce res

time = mea.Time();
S = mea.seizure;
onset_inds = find(diff(S(:, 3)) > 0) + 1;
onset_times = time(onset_inds);
nearinds = interp1(onset_times, 1:length(onset_times), res.time, ...
	'nearest', 'extrap');

onset_times = onset_times(nearinds);
onset_inds = onset_inds(nearinds);

gt.computeTimes = res.time * 1e3;
gt.Z = angle(-S(onset_inds, 1:2) * [1; 1j]);
gt.speed = abs(S(onset_inds, 1:2) * [1; 1j]) ./ (res.time - onset_times');
gt.V = [gt.speed .* cos(gt.Z), gt.speed .* sin(gt.Z)]';
gt.p = 0 * nearinds;
gt.beta = pinv(gt.V');

fits.groundtruth = gt;

end