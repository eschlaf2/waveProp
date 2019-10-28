function [new_map,state] = update_map(params, state)

Nx = params.grid_size(1);
Ny = params.grid_size(2);
map_type = params.map_type;

%center of initial source map.
xCenter = params.stim_center(1); 
yCenter = params.stim_center(2); 

if ~exist('state', 'var'), state = []; end

if strcmp(map_type, 'fixed_point_source')
	new_map = false(Nx,Ny);

	%set the initial map.
	new_map(xCenter + (-1:1), yCenter + (-2:0)) = 1;
	state = NaN;                                %Not used in this case.
	return
end

%% (map_type: 'ictal_wavefront')

shrink_factor = 0.5;
if numel(state) < Nx * Ny
	%set the initial map.
	state = false(Nx, Ny);
	state(xCenter + (-1:1), yCenter + (-1:1)) = 1;
	new_map = state; 
	return
end

[r,c] = find(state);
B = boundary(c,r,shrink_factor);

% For each point on the boundary,
for i=1:length(B)-1
	%Get a boundary point,
	rON = r(B(i));%B(i,1);
	cON = c(B(i));%B(i,2);
	%Get the neighbors,
	[xx, yy] = ndgrid(1:Nx, 1:Ny);
	isa_neighbor = ( (xx - rON).^2 + (yy - cON).^2 ) <= 2;
	candidates = (state == 0) & isa_neighbor;
	r_near = xx(candidates);
	c_near = yy(candidates);

	if ~isempty(r_near)
		ind0 = randi(length(r_near));
		state(r_near(ind0), c_near(ind0)) = 1;
	end
end

new_map = false(size(state));
B0 = bwboundaries(state, 'noholes');
B0 = B0{1};
for m=1:length(B0)
	new_map(B0(m,1),B0(m,2))=1;
end

end
