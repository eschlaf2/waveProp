function [map,state] = update_map(state, expansion_rate)

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

% shrink_factor = 0.5;
% [Nx, Ny] = size(state);
% [xx, yy] = ndgrid(1:Nx, 1:Ny);

% [r,c] = find(state);
% mask = state(:);
% D = squareform(pdist([xx(:), yy(:)]));
% neighbor_mask = any(D(mask, :) <= sqrt(2));

map = false(size(state));
map(2:end-1, 2:end-1) = conv2(state, [0 1 0; 1 -4 1; 0 1 0], 'valid') > 0;
if expansion_rate <= 0, return; end
% map = neighbors;
% neighbors(neighbor_mask) = true;
p_thresh = 2^expansion_rate - 1;
recruitment = ( (rand(size(state)) > (1 - p_thresh)) & map );
state = recruitment | state;
% map_decay = (rand(size(map)) > (1 - p_thresh)) & map;
% map = xor(map, map_decay) | xor(new_state, state);
% state = new_state;


% B = boundary(c,r,shrink_factor);
% 
% % For each point on the boundary,
% for i=1:length(B)-1
% 	
% 	%Get a boundary point,
% 	rON = r(B(i));%B(i,1);
% 	cON = c(B(i));%B(i,2);
% 	
% 	%Get the non-recruited neighbors,
% 	isa_neighbor = ( (xx - rON).^2 + (yy - cON).^2 ) <= 2;  % check proximity
% 	candidates = (state == 0) & isa_neighbor;  % ... and state
% 	r_near = xx(candidates);
% 	c_near = yy(candidates);
% 
% 	% Choose a random neighbor and activate it
% 	if ~isempty(r_near)
% 		ind0 = randi(length(r_near));
% 		state(r_near(ind0), c_near(ind0)) = 1;
% 	end
% end
% 
% Set the map to be the boundary between recruited and non-recruited territory 
% map = false(size(state));
% B0 = bwboundaries(state, 'noholes');
% B0 = B0{1};
% map(sub2ind(size(map), B0(:, 1), B0(:, 2))) = true;


end
