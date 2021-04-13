function [f, xi, bw] = circ_ksdens(theta, varargin)
     % Wrapper to compute the mode of angular data using ksdensity
    % Input <theta> should be given in radians
    % Default bandwidth is pi/12; default resolution is pi/64;
    % takes arguments after <theta> are passed to ksdensity

    BW = pi/12;
    XI = linspace(-pi, pi, 129);
    
    if isvector(theta), theta = theta(:); end
    
    data = [theta + 2*pi; theta; theta - 2*pi];
    
    if ~isempty(varargin) && isnumeric(varargin{1})  % if xi is given in varargin ...
        XI = varargin{1};
        varargin(1) = [];
        [f, xi, bw] = ksdensity(data, XI, 'bandwidth', BW, varargin{:});
    else
        [f, xi, bw] = ksdensity(data, XI, 'bandwidth', BW, varargin{:});
    end
    
end
