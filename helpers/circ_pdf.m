function [hist, centers] = circ_pdf(dirs, centers, range, varargin)
% Smoothed histogram of circular data over <range>

    defaultCenters = linspace(-pi, pi, 129);
    if nargin < 2 || isempty(centers), centers = defaultCenters; end
    if ischar(centers), 
        varargin = [centers range varargin]; 
        range = [];
        centers = defaultCenters;
    end
    if nargin < 3 || isempty(range), range = [-pi pi]; end
    if ischar(range), varargin = [range, varargin]; end
    
    % Default name-value pair arguments
    P = inputParser;
    p = @(varargin) addOptional(P, varargin{:});
    
    p('PeakToTheta', [], @(x) isempty(x) || x >= range(1) && x <= range(2));  % shift highest peak to appear at angle PeakToTheta
    p('SmoothingMethod', 'gaussian', @(x) smoothdata(1, x));
    p('SmoothingWindow', pi/4, @(x) x >= 0 && x <= 2*pi);
    p('PeakSpecs', {}, @(x) iscell(x));
    p('Normalization', 'pdf');

    P.parse(varargin{:});
    p = P.Results;
    
    edges = centers - diff(centers(1:2))/2;
    counts = histcounts(dirs, edges, 'normalization', p.Normalization);
    
    centers = centers(1:end-1);
    
    cc = [centers-2*pi centers centers + 2*pi];
    counts_temp = [counts counts counts];
    counts_sm = smoothdata(counts_temp, ...
        p.SmoothingMethod, p.SmoothingWindow, 'SamplePoints', cc);
    
    % mask to original interval
    mask = abs(cc) <= pi;
    
    
    if ~isempty(p.PeakToTheta)
        
        % find peaks
        directives = ['NPeaks', 2, 'SortStr', 'descend', p.PeakSpecs];
        [~, locs] = findpeaks(counts_sm(mask), cc(mask), directives{:});
    
        cc = cc - (locs(1) - p.PeakToTheta);
    end
    
    hist = interp1(cc, counts_sm, centers, 'nearest');

    
end



