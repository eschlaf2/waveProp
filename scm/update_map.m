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
	state = map; 
	return
end

%% Update

if expansion_rate <= 0, map = state; return; end

map = false(size(state));
map(2:end-1, 2:end-1) = conv2(state, [0 1 0; 1 -4 1; 0 1 0], 'valid') > 0;
p_thresh = 2^expansion_rate - 1;
if exist('excitability_map', 'var')
	recruitment = ( (rand(size(state)) .* excitability_map > (1 - p_thresh)) & map );
	map = map & excitability_map;
else
	recruitment = ( (rand(size(state)) > (1 - p_thresh)) & map );
end
state = recruitment | state;


end
