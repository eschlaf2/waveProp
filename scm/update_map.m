function [map,state] = update_map(state, expansion_rate, excitability_map, ~)

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

if expansion_rate <= 0, map = logical(state); return; end

% boundary = zeros(size(state));
boundary = conv2(state > 0, [0 1 0; 1 -4 1; 0 1 0], 'same') > 0;
if ~exist('excitability_map', 'var'); excitability_map = ones(size(state)); end

p_wavefront = 2^(.5*expansion_rate) - 1;  % 50% of area is recruited at this rate
p_recruit = 2.^(.5*expansion_rate * excitability_map) - 1;  % wavefront dies off at rate indicated by excitability map
boundary = boundary .* excitability_map;

dice = rand(size(state));
wavefront = (dice < p_wavefront) & boundary & (state == 0);
recruit = dice < p_recruit & (state == 1);

state(wavefront) = 1;
state(recruit) = 2;
map = state == 1;


end
