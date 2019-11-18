function [map,state] = update_map(state, expansion_rate, excitability_map)

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
if exist('excitability_map', 'var')
	p_thresh = 2.^(expansion_rate * excitability_map) - 1;
	boundary = boundary .* excitability_map;
else
	p_thresh = 2^expansion_rate - 1;
end
dice = rand(size(state));
new_recruits = (dice < 2^expansion_rate - 1) & boundary & (state == 0);
ictal = dice < p_thresh & (state == 1);

state(new_recruits) = 1;
state(ictal) = 2;
map = state == 1;


end
