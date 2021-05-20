classdef MaxDescent < WaveProp
	
	properties
		HalfWin = 0.05  % seconds
		FBand = [1 50]  % Hz
        Ascent = false  % use max *ascent* (instead of max *descent*)
        MinFiniteLow = 10
        DiffsOrPeaks string {mustBeMember(DiffsOrPeaks,{'diffs','peaks'})} = 'diffs'
    end
    
    properties (Transient = true)
        % Used to get the TOA, but no need to save afterward
        Time  % This will just be indices unless given as input
        SamplingRate = 1  % Assumes 1 sample per second unless otherwise given
    end
	
	
	methods
        
		function obj = MaxDescent(signal, t0, varargin)
			% Get the TOA at times <t0> in the given <signal>. If you give
			% a signal directly, 'Time' and 'SamplingRate' will be in
			% indices so you will need to give those as arguments to change
			% 'HalfWin'
            
            % Parse the inputs
			if nargin < 1, return, end  % return an empty object
			if ischar(signal), obj = obj.parse_inputs(signal, t0, varargin{:}); return; end
            if isa(signal, 'MEA')  % parse from MEA
				mea = signal;
				obj.t0 = t0(:);
				obj.HalfWin = mea.params.half_win / 1e3;
				obj.Position = mea.Position;
                obj.Name = mea.Name;
                obj.Time = mea.Time;
                obj.GridSize = mea.GridSize;
                obj.SamplingRate = mea.SamplingRate;
				obj = obj.parse_inputs(varargin{:});
                if isempty(mea.MaxDescentData)
					mea.MaxDescentData = zscore(mea.filter(mea.Data, mea.SamplingRate, obj.FBand));
                end
                signal = mea.MaxDescentData;
            else  % Return the TOA based on the entire signal
                obj.t0 = t0(:);
				obj = obj.parse_inputs(varargin{:});
                if ndims(signal) == 3  % Get positions and reshape to 2D
                    obj.GridSize = size(signal, [2 3]);
                    [px, py] = ndgrid(1:obj.GridSize(1), 1:obj.GridSize(2));
                    obj.Position = [px(:) py(:)];
                    signal = signal(:, :);
                end
            end
			
            if isempty(obj.Time), obj.Time = 1:size(signal, 1); end
            
            
            % flip the sign of the signal if looking for Ascents
            if obj.Ascent, signal = -signal; end
            
            % Get the TOA
			obj.TOA = obj.get_toa(signal);
            if obj.HalfWin < .05, obj.MinFinite = obj.MinFiniteLow; end
			
            % Get rid of this... (i.e. do all computation in WaveProp)
%             obj = obj.compile_results;
			
        end
        
        function toa = get_toa(obj, signal)
            % Find descent/peak times within HalfWin of t0. Result is in
            % time from window start (i.e. t0 - HalfWin)
            
            % Make sure position matches signal
            assert(size(obj.Position, 1) == size(signal, 2));
            
            % Initialize toa
            toa = nan(numel(obj.t0), size(signal, 2));
            
            % Get the toa for the window surrounding each t0
            for ii = 1:numel(obj.t0)
                % Get the window surrounding t0
                inds = abs(obj.Time - obj.t0(ii)) <= obj.HalfWin;
                window = signal(inds, :);
                
                window = window - window(1, :);  % set first time point as baseline (for visualization early)
			
                if strcmpi(obj.DiffsOrPeaks, 'peaks')
                    [change, time_point] = min(window);  % Find time of minimal peak
                else
                    [change, time_point] = min(diff(window, 1, 1));  % Find time of maximal descent
                end
                non_decreasing = change >= 0;  % Find non-decreasing traces
                bdry = (time_point == 1) | (time_point >= size(window, 1) - 1);  % ... and traces with max descent on the boundary (these are not part of the wave and confuse the analysis)
                inactive = range(window) < 1;

                data = time_point;
                data(non_decreasing | bdry | inactive) = nan;
                toa(ii, :) = data ./ obj.SamplingRate;
            end
            
            
        end
        
        
        function obj = reload(obj, S)
            obj = reload@WaveProp(obj, S);
        end
    end       

    methods (Static)
        function obj = loadobj(S)
            if isstruct(S)
                obj = MaxDescent;
                for ff = string(fieldnames(S)')
                    if ismember(ff, {'ParamNames', 'WPParamNames'}), continue; end
                    obj.(ff) = S.(ff); 
                end
                obj = reload(obj, S);
            else
                obj = S;
            end
        end
    end
        
%% Old 
% 		function [window, t] = get_window(obj, mea)
% 			t_inds = abs(mea.Time - obj.t0) <= obj.HalfWin;
% 			t = mea.Time(t_inds);
% 			data = mea.MaxDescentData;
%             if obj.Ascent, data = -data; end
% 			window = data(t_inds, :);
%         end
%         
% 		function data = get_data(M, window)
% 			% Get the TOAs from the given window
%             
% 			window = window - window(1, :);  % set first time point as baseline (for visualization early)
% 			
%             if strcmpi(M.DiffsOrPeaks, 'peaks')
%                 [change, time_point] = min(window);  % Find time of minimal peak
%             else
%                 [change, time_point] = min(diff(window, 1, 1));  % Find time of maximal descent
%             end
% 			non_decreasing = change >= 0;  % Find non-decreasing traces
% 			bdry = (time_point == 1) | (time_point == size(window, 1) - 1);  % ... and traces with max descent on the boundary (these are not part of the wave and confuse the analysis)
% 			inactive = range(window) < 1;
% 			
% 			data = time_point;
% 			data(non_decreasing | bdry | inactive) = nan;
% 			
%         end
%         
% 	end
	
	
	
end