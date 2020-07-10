classdef WaveFlow < WaveProp
	
	properties
		flow
		data_s
	end
	
	properties (Access = public)  % set to private
		K = 3
		smoothing_window
		scale_factor = 5
		sampling_rate
		min_flow_speed = 40  % in mm/s ... (cross the MEA in < .1s)
	end
	
	properties (Dependent = false)  % set to true so you don't retain these later
		rescale
% 		data_i
% 		data_o
% 		Zc
% 		Vx
% 		Vy
% 		Direction
% 		p
	end
	
	methods
		function obj = WaveFlow(mea, times)
			if nargin < 2, times = []; end
			if nargin < 1, mea = MEA; end
			obj.sig = -obj.min_flow_speed;
% 			obj.data_o.data = mea.Data;
% 			obj.data_o.position = mea.Position;
% 			obj.sampling_rate = mea.SamplingRate;
			obj.smoothing_window = [3, 3, max(1 + round(10 * mea.SamplingRate / 1e3), 3)];
% 			obj.Time = mea.Time;
			[window, obj.t0] = obj.get_window(mea.Data, mea.Time, times);
			data_i = WaveFlow.interpolate_nans(window, mea.Position);
			obj.data_s = WaveFlow.smooth_data(data_i, obj.smoothing_window, obj.K);
			obj.flow = WaveFlow.compute_flow(obj.data_s);
			obj.flow.Magnitude = obj.rescale * obj.flow.Magnitude;
			obj.flow = WaveFlow.mask_data(obj.flow, obj.K, obj.min_flow_speed);
			obj = obj.summarize;
		end
		
		function s = summarize(s)
			vx = s.flow.Vx;
			vy = s.flow.Vy;
			s.Vx = squeeze(nansum(vx, [1 2]));
			s.Vy = squeeze(nansum(vy, [1 2]));
			N = squeeze(sum(isfinite(complex(vx, vy)), [1 2]));
			Z = squeeze(nansum(complex(vx, vy), [1 2]));
% 			Z = reshape(complex(vx, vy), [], size(vx, 3));
			
			s.Magnitude = s.rescale * abs(Z) ./ N;
			dir = angle(Z);
			dir(isnan(s.Magnitude)) = nan;
			s.Direction = dir;
			s.p = -s.Magnitude;
		end
		function r = get.rescale(s)
			r = .4 / 2.^s.K / diff(s.t0(1:2));
		end
		
		function ax = plot(s, varargin)
			ax = []; t0 = s.t0(1);
			for arg = varargin
				switch class(arg{:})
					case 'matlab.graphics.axis.Axes'
						ax = arg{:};
				
					case {'double', 'single'}
						t0 = arg{:};
					otherwise
						error('Argument not recognized. Must be numeric or Axes.')
				end
			end
			if isempty(ax), ax = gca; end
			[d, which_t] = min(abs(s.t0 - t0));
			if d > 1e-3, return; end
			imagesc(ax, s.data_s(:, :, which_t));
			colormap(ax, 1-gray);
% 			colormap bone;
			hold(ax, 'on')
			quiver(ax, s.flow.Vx(:, :, which_t), s.flow.Vy(:, :, which_t), s.scale_factor);
			hold(ax, 'off')
		end
		

	end
	
	methods (Static)
		function [window, t0] = get_window(data, time, times)
			if isempty(times)
				window = data;
				t0 = time;
			elseif numel(times) == 2
				t_inds = time > time_range(1) & time < time_range(2);
				window = data(t_inds, :);
				t0 = time(t_inds);
			else
				t_inds = interp1(time, 1:length(time), times, 'nearest', 'nearest');
				window = data(t_inds, :);
				t0 = time(t_inds);
			end
		end
		
		function data_i = interpolate_nans(data, position)
			disp('Interpolating ...')
			data = -zscore(single(data));
			Nt = length(data);
			layout = max(position);
			Nr = layout(1); Nc = layout(2);
			data_i = nan(Nr, Nc, Nt);
% 			[xx, yy] = ndgrid(linspace(1, Nr, Nr), linspace(1, Nc, Nc));
			[xx, yy] = ndgrid(1:Nr, 1:Nc);
			for ii = 1:Nt
				F = scatteredInterpolant(position(:, 1), position(:, 2), ...
					double(data(ii, :))', 'natural'); 
				data_i(:, :, ii) = reshape(F(xx(:), yy(:)), Nr, Nc); 
			end
		end
		function data_s = smooth_data(data, smoothing_window, K)
			disp('Smoothing ...')
% 			dataS = del2(smooth3(data, 'gaussian', smoothing_window));
			data = smooth3(data, 'gaussian', smoothing_window);
			nT = size(data, 3);
% 			[Xq, Yq] = meshgrid(1:size(data, 1), 1:size(data, 2));
% 			Xq = interp2(Xq, K);
% 			Yq = interp2(Yq, K);
			for ii = nT:-1:1
				temp = del2(data(:, :, ii));
				data_s(:, :, ii) = ...
					single(interp2(temp, K, 'spline'));
			end
			data_s = data_s - mean(data_s, 3);
		end
		function flow = compute_flow(data_s)
			fprintf('Computing flow ... ')
			opticFlow = opticalFlowLK;
			nT = size(data_s, 3);
			flow(nT) = opticalFlow;
			for ii = 1:nT
				flow(ii) = estimateFlow(opticFlow, data_s(:, :, ii));
			end
			flow = WaveProp.resize_obj(flow, 1);
			fprintf('Done.\n')
			
		end
		
		function flow = mask_data(flow, K, min_wave_speed)
			too_slow = flow.Magnitude < min_wave_speed;
			bdy = true(size(flow.Vx));
			margin = 2^(K-1);
			bdy(margin+1:end-margin, margin+1:end-margin, :) = false;
			mask = too_slow | bdy;
			for f = string(fieldnames(flow)')
				flow.(f)(mask) = nan;
			end
		end

	end
	
end