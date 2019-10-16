% default_params
function [res, G] = init_params(varargin)

%% parsing functions
G = inputParser;  % Define the parser
p = @(varargin) addParameter(G, varargin{:});  % convenience function
validate = @(x, all) any(validatestring(x, all));  % validate strings
isnonneg = @(x) isnumeric(x) & x >= 0;  % check nonnegative

%% wave_prop parameters
allFitMethods = {'nyc', 'bos'};

p('fit_method', 'bos', @(x) validate(x, allFitMethods));
p('T', 10);  % window length for delay method (s)
p('half_win', 50); % half-window for all non-delay methods (ms)
p('exclude', false, @islogical);  % logical: exclude inactive channels?
p('thresh', Inf);  % threshold for deviance methods
p('delay_band', [1 13]);  % frequency band for coherence
p('show_plots', false, @islogical);  % show detailed plots during analysis



%% mua_events
p('event_thresh', 4, isnonneg);      % min peak distance from baseline (in negative direction) [sd]
p('min_dist', 1, isnonneg);          % min time between peaks [ms]
p('artefact_thresh', 16, isnonneg);  % activity outside this range is considered artefact [sd] 

%% get_discharge_times
p('discharge_method', 1, @(x) any(x == [1, 2]));  % discharge times based on (1) peaks in firing rate; (2) peaks in LFP

% Method 1 (firing rate based)
p('min_fr', 1/60);   % exclude channels with mean firing rate less than one spike per minute (Liou et al., 2018) ?\cite{Liou2018a}
p('min_peak_distance', 100);  % (ms)
p('min_peak_prominence', 100);  % (std of diff);

% Method 2  (lfp based)
p('min_lfp_deviance', 2);  % (sd)
p('min_active_electrodes', 10);  % (electrodes)
p('smoothing_win', 10);  % (ms)

%% exclude_channels
p('exclude_rate', 6);  % (spikes/sec) Exclude channels with mean firing rate less than this
p('require_increase', true, @islogical);  % Exclude channels that don't increase in firing rate during the first half of the seizure


%% Parse
parse(G, varargin{:})
res = G.Results;

end