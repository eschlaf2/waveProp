function [f, xi, bw] = circ_ksdens(theta, varargin)
     % Wrapper to compute the mode of angular data using ksdensity
    % Input <theta> should be given in radians
    % Default bandwidth is pi/12; default resolution is pi/64;
    % takes arguments after <theta> are passed to ksdensity

    BW = pi/12;
    XI = linspace(-pi, pi, 129);
    
    if isvector(theta), theta = theta(:); end
    
    data = [theta + 2*pi; theta; theta - 2*pi];
    [f, xi, bw] = ksdensity(data, XI, 'bandwidth', BW, varargin{:});
    
end
