function [map,state] = update_map_smooth(state, expansion_rate, excitability_map, dt)

%% Initialization
if isstruct(state) 
	Nx = state.grid_size(1); Ny = state.grid_size(2);

	%center of initial source map.
	xCenter = state.stim_center(1); 
	yCenter = state.stim_center(2); 

	map = false(Nx,Ny);

	%set the initial map.
% 	map(xCenter + (-1:1), yCenter + (-1:1)) = 1;
	map(xCenter, yCenter) = 1;
	state = zeros([Nx, Ny]);
	state(map) = 1;
	return
end

%% Update

if expansion_rate <= 0, map = logical(state >= 0); return; end

% boundary = zeros(size(state));
boundary = conv2(state > 0, [0 1 0; 1 -4 1; 0 1 0], 'same') > 0;
if ~exist('excitability_map', 'var'); excitability_map = ones(size(state)); end

p_wavefront = 2^(expansion_rate) - 1;  % area is recruited at this rate
% boundary = boundary .* excitability_map;

dice = rand(size(state));
wavefront = (dice < p_wavefront) & boundary & (state < 0);

recruited = state >= 0;
state(recruited) = state(recruited) + dt;

state(wavefront) = 0;
map = state >= 0 & state < excitability_map;


end
