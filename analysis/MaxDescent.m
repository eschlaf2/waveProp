classdef MaxDescent < WaveProp
	
	properties
		HalfWin = 0.05
		FBand = [1 50]
        Ascent = false
        MinFiniteLow = 10
        DiffsOrPeaks char {mustBeMember(DiffsOrPeaks,{'diffs','peaks'})} = 'diffs'
	end
	
	properties (Hidden = true)
		ParamNames = ["Position" "HalfWin" "FBand" "Ascent"]
	end
	
	methods
        
		function obj = MaxDescent(data, t0, varargin)
			
			if nargin < 1, return, end
			if ischar(data), obj = obj.parse_inputs(data, t0, varargin{:}); return; end
            if isa(data, 'MEA')
				mea = data;
				obj.t0 = t0;
				obj.HalfWin = mea.params.half_win / 1e3;
				obj.Position = mea.Position;
                obj.Name = mea.Name;
				obj = obj.parse_inputs(varargin{:});
				if isempty(mea.MaxDescentData)
					mea.MaxDescentData = zscore(mea.filter(mea.Data, mea.SamplingRate, obj.FBand));
				end
				window = obj.get_window(mea);
            else
				obj = obj.parse_inputs(obj.ParamNames, varargin{:});
				window = data;
            end
			
	
			obj.Data = obj.get_data(window);
			obj.Data = obj.Data / mea.SamplingRate;
            if obj.HalfWin < .05, obj.MinFinite = obj.MinFiniteLow; end
			obj = obj.compile_results;
			
        end
        
		function [window, t] = get_window(obj, mea)
			t_inds = abs(mea.Time - obj.t0) <= obj.HalfWin;
			t = mea.Time(t_inds);
			data = mea.MaxDescentData;
            if obj.Ascent, data = -data; end
			window = data(t_inds, :);
        end
        
		function data = get_data(M, window)
			
			window = window - window(1, :);  % set first time point as baseline (for visualization early)
			
            if strcmpi(M.DiffsOrPeaks, 'peaks')
                [change, time_point] = min(window);  % Find time of minimal peak
            else
                [change, time_point] = min(diff(window, 1, 1));  % Find time of maximal descent
            end
			non_decreasing = change >= 0;  % Find non-decreasing traces
			bdry = (time_point == 1) | (time_point == size(window, 1) - 1);  % ... and traces with max descent on the boundary (these are often not part of the wave and confuse the analysis)
			inactive = range(window) < 1;
			
			data = time_point;
			data(non_decreasing | bdry | inactive) = nan;
			
        end
        
	end
	
	
	
end