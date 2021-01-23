classdef Events < WaveProp
	
	properties
		HalfWin = .05
	end
	
	properties (Hidden = true)
		ParamNames = "HalfWin"
	end
	
	methods
		
		function obj = Events(mea, t0, varargin)
			
			if nargin < 1, return, end
			assert(isa(mea, 'MEA'));
% 			obj.Time = mea.AllTime;
			obj.t0 = t0;
			obj.HalfWin = mea.params.half_win / 1e3;
			obj = obj.parse_inputs(varargin{:});
			obj.Position = mea.Position;
            obj.Name = mea.Name;
% 			[window, T] = obj.get_window(mea);

			data = obj.get_data(mea.event_times, t0);
			obj.Data = data;
			obj = obj.compile_results; 			
		end
		function [window, T] = get_window(obj, mea)  % Not used
			t_inds = abs(mea.AllTime - obj.t0) <= obj.HalfWin;
			events = mea.mua_events(t_inds, :);
			window = events;
			T = mea.AllTime(t_inds);
		end
		function data = get_data(obj, event_times, t0)
			
			[ch, ET] = event_times{:};
			mask = abs(ET - t0) <= obj.HalfWin;
			ch(~mask) = [];
			ET(~mask) = [];
			N = size(obj.Position, 1);
			data = nan(N, 1);
			for ii = 1:N
% 				ET = event_times{ii};
% 				ET(abs(ET - t0) > obj.HalfWin) = [];
				data(ii) = nanmean(ET(ch == ii));
			end
			
		end
	end
	
	
	
end