function [dens, xi] = circ_ksdens(dirs, xi, bw)
    % [dens, xi] = circ_ksdens(dirs, xi=(-pi:pi/128:pi), bw=pi/16)
    % Circular ksdensity. <xi> and <bw> are used as in ksdensity().
    % Inputs: <dirs> will be converted to Nx1.
    %    dirs: Nx1 vector of directions in radians
    %    xi: resolution (scalar) or xi (vector; see <ksdensity.m>)
    %    bw: bandwidth (scalar) (see <ksdensity.m>)
    
    if nargin < 2 || isempty(xi), xi = 128; end
    if nargin < 3, bw = pi/16; end
    
    dirs = dirs(:);
    if numel(xi) == 1
%         xi = linspace(-pi/2, 3*pi/2, xi);
        xi = linspace(-pi, pi, xi);
    end
    
    d0 = dirs + [2*pi 0 -2*pi];
    dens = ksdensity(d0(:), xi, 'bandwidth', bw);
end
