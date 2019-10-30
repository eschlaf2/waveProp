function [map,state] = update_map(state)

%% Initialization
if isstruct(state) 
	Nx = state.grid_size(1); Ny = state.grid_size(2);

	%center of initial source map.
	xCenter = state.stim_center(1); 
	yCenter = state.stim_center(2); 

	map = false(Nx,Ny);

	%set the initial map.
	map(xCenter + (-1:1), yCenter + (-1:1)) = 1;
	state = map;                                %Not used in this case.
	return
end

%% Update

shrink_factor = 0.5;
[Nx, Ny] = size(state);
[xx, yy] = ndgrid(1:Nx, 1:Ny);

[r,c] = find(state);
B = boundary(c,r,shrink_factor);

% For each point on the boundary,
for i=1:length(B)-1
	
	%Get a boundary point,
	rON = r(B(i));%B(i,1);
	cON = c(B(i));%B(i,2);
	
	%Get the non-recruited neighbors,
	isa_neighbor = ( (xx - rON).^2 + (yy - cON).^2 ) <= 2;  % check proximity
	candidates = (state == 0) & isa_neighbor;  % ... and state
	r_near = xx(candidates);
	c_near = yy(candidates);

	% Choose a random neighbor and activate it
	if ~isempty(r_near)
		ind0 = randi(length(r_near));
		state(r_near(ind0), c_near(ind0)) = 1;
	end
end

% Set the map to be the boundary between recruited and non-recruited territory 
map = false(size(state));
B0 = bwboundaries(state, 'noholes');
B0 = B0{1};

for m=1:length(B0)
	map(B0(m,1), B0(m,2))= 1;
end

end
