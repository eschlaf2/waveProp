function [mode, f, xi, bw] = circ_mode(theta, varargin)
    % Wrapper to compute the mode of angular data using ksdensity
    % Input <theta> should be given in radians
    % Default bandwidth is pi/12; default resolution is pi/64;
    % takes arguments after <theta> are passed to ksdensity
    if isempty(theta); mode = nan; return; end
    BW = pi/12;
    XI = linspace(-pi, pi, 129);
    
    if isvector(theta), theta = theta(:); end
    
    data = [theta + 2*pi; theta; theta - 2*pi];
    [f, xi, bw] = ksdensity(data, XI, 'bandwidth', BW, varargin{:});
    
    [~, mode_i] = max(f);
    mode = xi(mode_i);
    
end