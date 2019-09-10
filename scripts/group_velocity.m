% group_velocity
% basename = compute_coherograms(pat, seizure, T, W, DS, units, toi);  %
% (See quickscript)
load(basename, ...  % Load the following from basename
    'f', 'params', 'phi', 'C', 't', 'confC', 'pat', ...
    'seizure', 'position', 'pairs', 'units')

df = mean(diff(f)); 
winsz = 3 * units;  % Hz
thresh = 5e-2; 
fband = params.fpass;
if ~exist('toi', 'var'); toi = [-Inf Inf]; end
if ~exist('tau', 'var'); tau = 'group'; end  % 'group' or 'phase'
% MASK = false;

if isinteger(phi); phi = -single(phi) / 1e4; end
finds = (f >= fband(1)) & (f <= fband(2));
tinds = (t >= toi(1)) & (t <= toi(2));
t = t(tinds); f = f(finds);

[nt, nf, np] = size(C(tinds, finds, :));

transform = @(A) reshape(permute(A(tinds, finds, :), [2 1 3]), nf, []);
Cf = transform(C);  % limit to band of interest
phif = transform(phi);  % ... same for phi

