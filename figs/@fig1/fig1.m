
classdef fig1 < BVNY & handle
	
	methods 
        function args = ksargs(F)
            
            args.Phi_mean_early = {linspace(-pi, pi, 100)};
            args.Phi_mean_late = {linspace(-pi/10, pi/10, 100)};
            args.Phi_std_early = {linspace(0, pi/2, 100), 'Support', [0 pi], 'boundarycorrection', 'reflection'};
            args.Phi_std_late = {linspace(0, pi/16, 100), 'Support', [0 pi], 'boundarycorrection', 'reflection', 'bandwidth', pi/512};
            args.First_detection = {linspace(0, 30, 100), 'support', [0 60], 'boundarycorrection', 'reflection'};
            args.N_detections_early = {'support', [-10 1e5], 'boundarycorrection', 'reflection'};
            args.Phi_early = {linspace(-pi, pi, 100)};
            
            if strcmpi(F.SimType, 'fs')
                args.Phi_late = {linspace(-pi/4, pi/4, 100), 'bandwidth', pi/64};
            else
                args.Phi_late = {'bandwidth', pi/64};
            end

        end
	end
	methods
		
        function F = fig1(varargin)
            F.AllSubs = 'ABDEG';  % C was events method; H was D1 method; F was first detection time
			if nargin < 1, return, end
            F.make;
            F.save;
            F.close;
        end
        
        save_fr_dir(F, sub)
		fr_dir_plot(F, sub)
		load_fits(F)
        [ax, data] = stats_summary(F, summary_metric)
        t = check_directions(N, metric)
        
        function h = panel_fig(F)
            fname = sprintf('SCM/%s/%s_%d_info.mat', F.SimType, F.SimType, F.DetailSim);
            load(fname, 'params', 'Qe_movie');
            h = figure(); fullwidth;
            set(h, 'units', 'points', 'position', [30   500   939   294])
            set(h, 'defaultAxesFontSize', 13)

            t = tiledlayout(h, 3, 10);
            set(t, 'padding', 'compact', 'tilespacing', 'compact');

            ctr = params.electrodes.centerNP;
            ext = params.electrodes.dimsNP;
            % stim = params.model.stim_center;
            stim = nan(size(params.source, 3), 2);
            for ii = 1:size(params.source, 3)
                [xi, xj] = find(params.source(:, :, ii) == max(params.source(:, :, ii), [], 'all'));
                stim(ii, :) = mean([xi xj]);
            end
            x = ctr(1) + [1 ext(1) ext(1) 1] - round(ext(1)/2);
            y = ctr(2) + [1 1 ext(2) ext(2)] - round(ext(2)/2);
            np =@(ax) patch(ax, x, y, 'red', 'linestyle', 'none', 'facealpha', .5);
            for ii = 11:40
                ax = nexttile(t, ii-10); 
                contourf(ax, Qe_movie(ii).cdata, 'linestyle', 'none'); 
                colormap(ax, 1 - gray); 
                np(ax); 
                hold(ax, 'on')
                plot(ax, stim(:, 2), stim(:, 1), 'r.', 'markersize', 15);
                hold(ax, 'off')
                title(num2str(ii-10)); 
            end
            axis(t.Children, 'square');
            set(t.Children, 'xtick', [], 'ytick', []);
            % set(h, 'position', [0    0.5633    0.6521    0.3222])
            h.Tag = 'scm_panel';
            F.A = h;
        end
        
        function args = KSArgs(F)
            args = F.ksargs();
        end
        
        function phi = true_angle(F, ii)
            fname = sprintf('SCM/%s/%s_%d_info.mat', F.SimType, F.SimType, F.WhichSims(ii));
            load(fname, 'params');
            phi = params.theta_FS;
            phi = phi(1);
%             phi = reshape(params.theta_FS, 1, []);
            
        end
        
		function save(F, subs)
            if nargin < 2 || isempty(subs), subs = F.AllSubs; end
            for sub = upper(subs)
                if isempty(F.(sub)), continue; end
               switch sub
                   case 'A'
                       save_A_(F);
                   case {'B', 'C', 'D', 'H'}
                       save_fr_dir(F, sub);
                   case 'E'
                       if strcmpi(F.SimType, 'sw')
                           F.save_EFG(F.E, [-90 180]);
                       else
                           F.save_EFG(F.E);
                       end
                   case 'F'
                       F.save_EFG(F.F);
                   case {'G'}
%                        temp = F.(sub);
                       for ii = 1:2
                        if F.(sub)(ii) == 0, continue; end
                        F.save_EFG(F.(sub)(ii));
                                
                       end
               end
           end
        end
			
			
        function close(F, subs)
            if nargin < 2 || isempty(subs), subs = F.AllSubs; end
            for sub = upper(subs)
                if isempty(F.(sub)), continue, end
                close(F.(sub));
                F.(sub) = [];
            end
        end
        function make(F, subs)
            if nargin < 2 || isempty(subs), subs = F.AllSubs; end
            for sub = upper(subs)
               switch sub
                   case 'A'
                       F.panel_fig;
                       
                   case {'B', 'D'}  % C and H were E and D10* methods
                       F.fr_dir_plot(sub);
                       
                   case 'ZE'  
                       ax = F.stats_summary("Phi_early");
                       title(ax, sprintf('Early (t < %d)', F.Early(2)));
                       F.E(2) = ax.Parent;
                       
                       ax = F.stats_summary("Phi_late");
                       title(ax, sprintf('Late (t > %d)', F.Late(1)));
                       F.E(1) = ax.Parent;
                       
                   case 'E'
                       F.Late(1) = -Inf;  % No distinction between early and late anymore
                       ax = F.stats_summary("Phi_late");
                       F.E(1) = ax.Parent;
                       
                   case 'ZF'
                       ax = F.stats_summary("First_detection");
                       title(ax, 'First Detection Time');
                       F.F = ax.Parent;
                       
                   case 'G'
                       ax = F.stats_summary("Phi_std_late");
                       title(ax, sprintf('Late (t > %d)', F.Late(1)));
                       F.G(2) = ax.Parent;
                       
                   case 'Z'
                       ax = F.stats_summary("Phi_std_early");
                       title(ax, sprintf('Early (t < %d)', F.Early(2)));
                       F.G(1) = ax.Parent;
               end
            end
        end

	end
	methods  % getters
		function files = get.Files(F)
            if isempty(F.Files)
                fname = @(ii) sprintf('%s_Seizure%d', F.SimType, ii);
                F.Files = arrayfun(fname, F.WhichSims, 'uni', 0);
            end
            files = F.Files;
        end
        function fits = get.Fits(F)
           if isempty(F.Fits), F.load_fits; end
           fits = F.Fits;
        end
        function pfx = prefix(F)
            if strcmpi(F.SimType, 'fs')
                pfx = [F.Prefix 'fig1/fig1_'];
            else
                pfx = [F.Prefix 'fig2/fig2_'];
            end
        end
        
        
    end
    
	properties
% 		Mea = MEA('FS/FS_Seizure101_Neuroport_10_10.mat');
		WhichSims = 1:10  %101:200
		SimType = 'FS'
        DetailSim = 9
		Early = [-Inf, 25]
		Late = [30 Inf]
		Files
		Fits
		Metrics = struct(...
			'B', 'M', ... 
			'C', '', ...
			'D', 'D10', ...
			'H', '');
		
		A
		B
		C
		D
		E
		F
        G
        H
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%
%% Function definitions

function save_A_(F)
    savefig(F.A, [F.prefix 'a_scm_panel'])
    print(F.A, [F.prefix 'a_scm_panel'], '-dpng')
end











%%
% style = set_style_();
% % runsplit2fit();
% 
% % A
% % panel_fig_()
% 
% % BCDH
% % create_bcdh_(style);
% 
% % % EFG
% fits = load_efg_();
% ax = create_efg_(fits, style);
% % export_efg_(ax);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %% Function definitions
% 
% function t = panel_fig_()
% % load('SCM_302_info.mat')
% load('SCM/FS/FS_112_info.mat', 'params', 'Qe_movie');
% h = figure(); fullwidth;
% set(h, 'units', 'points', 'position', [30   500   939   294])
% set(h, 'defaultAxesFontName', 'helvetica', ...
% 	'defaultAxesFontSize', 13)
% 
% t = tiledlayout(h, 3, 10);
% set(t, 'padding', 'compact', 'tilespacing', 'compact');
% 
% ctr = params.electrodes.centerNP;
% ext = params.electrodes.dimsNP;
% stim = params.model.stim_center;
% x = ctr(1) + [1 ext(1) ext(1) 1] - round(ext(1)/2);
% y = ctr(2) + [1 1 ext(2) ext(2)] - round(ext(2)/2);
% np =@(ax) patch(ax, x, y, 'red', 'linestyle', 'none', 'facealpha', .5);
% for ii = 11:40
% 	ax = nexttile(t, ii-10); 
% 	contourf(ax, Qe_movie(ii).cdata, 'linestyle', 'none'); 
% 	colormap(ax, 1 - gray); 
% 	np(ax); 
% 	hold(ax, 'on')
% 	plot(ax, stim(2), stim(1), 'r*');
% 	hold(ax, 'off')
% 	title(num2str(ii)); 
% end
% axis(t.Children, 'square');
% set(t.Children, 'xtick', [], 'ytick', []);
% % set(h, 'position', [0    0.5633    0.6521    0.3222])
% savefig(h, 'figs/bos_v_ny/fig1/fig1_a_scm_panel')
% print(h, 'figs/bos_v_ny/fig1/fig1_a_scm_panel', '-dpng')
% 
% 
% end
% 
% function create_bcdh_(style)
% % This is insane, but FS 101 direction is too perfect... waves travel perfectly along a straight line and Events can't fit anything...
% metrics = {'M', 'E', 'D10', 'D1xwh'};
% mea = MEA('FS/FS_Seizure112_Neuroport_10_10.mat');
% % mea = MEA('SCM/SCM_Seizure302_Neuroport_10_10.mat');
% fits = WaveProp.load(mea, 'metrics', metrics);
% s = load('figs/bos_v_ny/fig1/fig1_efg_style');
% s.Width = 2.5;
% s.Height = 1;
% s.FixedLineWidth = '1';
% s.FixedFontSize = 9;
% % s.ApplyStyle = true;
% close all
% 
% for mm = string(metrics)
% 	ax.(mm) = axes(figure);
% 	plot(ax.(mm), mea.Time, rescale(mean(mea.firing_rate, 2), -pi, pi), 'color', .5*[1 1 1]);
% 	hold(ax.(mm), 'on');
% 	scatter(ax.(mm), fits.(mm).time, fits.(mm).Direction, 'filled', 'markerfacecolor', style.(mm).color, 'sizedata', 25);
% 	hold(ax.(mm), 'off');
% 	axis(ax.(mm), 'tight');
% 	ylim(ax.(mm), [-pi pi]);
% 	yticks(ax.(mm), [-pi 0 pi]);
% 	yticklabels(ax.(mm), {'-\pi', '0', '\pi'});
% 	h = ax.(mm).Parent;
% 	fname = ['figs/bos_v_ny/fig1/f1_bcdh_' mm{:}];
% 	savefig(h, fname); 
% 	
% % 	lgd = findobj(h, 'type', 'legend');
% % 	lgd.Visible = 'off';
% 	aa = ax.(mm);
% 	p0 = aa.Position;
% 	aa.Position = [0 0 1 1];
% 
% 	hgexport(h, fname, s);
% 	aa.Position = p0;
% end
% 
% end
% 
% function fits = load_efg_()
% which_sims = 101:160;
% % metrics = {'M', 'E', 'D10', 'D1', 'Mdt', 'Edt'};
% metrics = {'M', 'E', 'D10', 'D1xwh'};
% 
% N = length(which_sims);
% info_files(N) = dir(['SCM/FS/FS_' num2str(N) '_info.mat']);
% for ii = 1:N
% 	ss = which_sims(ii);
% 	info_files(ii) = dir(['SCM/FS/FS_' num2str(ss) '_info.mat']);
% end
% p_num = nan(size(info_files));
% for ii = 1:length(info_files)
% 	p_num(ii) = str2double(info_files(ii).name(4:end-9));
% end
% 
% clear files
% files(N) = dir(['FS_Seizure' num2str(ss) '_fits.mat']);
% f_num = which_sims;
% for ii = 1:N
% 	ss = which_sims(ii);
% 	files(ii) = dir(['FS_Seizure' num2str(ss) '_fits.mat']);
% % 	f_num(ii) = ss;
% end
% fits = WaveProp.load('files', files, 'metrics', metrics);
% 
% 
% true_angle = nan(N, 1);
% params =@(ii) load([info_files(ii).folder filesep info_files(ii).name], 'params').params;
% % [xx, yy] = ndgrid(1:50, 1:50);
% grid_size = params(1).model.grid_size;
% [xx, yy] = ind2sub(grid_size, 1:prod(grid_size));
% stim_center =@(params) [nanmean(xx(params.meta.source), 'all') nanmean(yy(params.meta.source), 'all')];
% get_angle =@(params) angle((params.electrodes.centerNP - stim_center(params)) * [1; 1j]);
% 
% 
% for ii = 1:N
% 	true_angle(ii) = get_angle(params(ii));
% end
% 
% for ii = 1:N
% 	for mm = string(metrics)
% % 		fits.(mm)(ii).Early = 20;
% % 		fits.(mm)(ii).Late = 30;
% 		fits.(mm)(ii).RotateBy = true_angle(ii);
% % 		if strcmpi(mm, 'D10')
% % 			fits.(mm)(ii).RotateBy = true_angle(ii) + pi;
% % 		else
% % 			fits.(mm)(ii).RotateBy = true_angle(ii);
% % 		end
% 	end
% end
% % fits.D1 = fits.D1xwh;
% % for ii = N:-1:1
% % 	fits.D10dt(ii) = fits.D10(ii).resample_t0(fits.Edt(ii).time); 
% % 	fits.D1dt(ii) = fits.D1(ii).resample_t0(fits.Edt(ii).time); 
% % end
% 
% 
% end
% 
% 
% function ax = create_efg_(fits, style)
% close all
% summ = struct;
% ax = struct;
% % ksargs.Phi_mean_early = {'numpoints', 3e3, 'bandwidth', .05};
% % ksargs.Phi_mean_early = {'numpoints', 150};
% ksargs.Phi_mean_early = {linspace(-pi, pi, 100)};
% % ksargs.Phi_mean_late = {'numpoints', 900, 'bandwidth', .25};
% % ksargs.Phi_mean_early = {linspace(-1.1*pi, 1.1*pi, 2000), 'bandwidth', .25};
% ksargs.Phi_mean_late = {linspace(-pi/10, pi/10, 100)};
% ksargs.Phi_std_early = {linspace(0, pi/2, 100), 'Support', [0 pi], 'boundarycorrection', 'reflection'};
% ksargs.Phi_std_late = {linspace(0, .1, 100), 'Support', [0 pi], 'boundarycorrection', 'reflection'};
% ksargs.First_detection = {linspace(0, 30, 100), 'support', [0 60], 'boundarycorrection', 'reflection'};
% ksargs.N_detections_early = {'support', [-1 1e5], 'boundarycorrection', 'reflection'};
% 
% 
% for ff = ["N_detections_early" "Phi_mean_early" "Phi_mean_late" "Phi_std_early" "Phi_std_late" "First_detection"]
% % for ff = ["N_detections_early" "Phi_mean_early" "Phi_mean_late"]
% 
% 	summ.(ff) = struct;
% 	ax.(ff) = axes(figure);
% % 	for mm = ["Mdt" "Edt" "D10dt" "D1dt" "M" "E" "D10" "D1"]
% 	for mm = string(fieldnames(fits)')
% 		if strcmpi(mm, 'name'), continue, end
% 		ww = strrep(mm, 'dt', '');
% 		col = style.(ww).color;
% 		if contains(mm, 'dt'), ls = '--'; else, ls = '-'; end
% 
% 		for ii = length(fits.Name):-1:1
% 			summ.(ff).(mm)(ii) = fits.(mm)(ii).(ff);
% 		end
% 		
% 		try
% 			data = summ.(ff).(mm);
% 			if contains(ff, 'Phi_mean')
% 				data = [data - 2*pi, data, data + 2*pi];
% 				data = data(abs(data) < pi * 1.25);
% 			end
% 			if isempty(data), fprintf('%s, %s bw: empty\n', ff, mm); continue; end
% 			[d, xi, bw] = ksdensity(data, ksargs.(ff){:});
% 			if contains(ff, 'Phi_mean')
% 				mask = xi > -pi & xi <= pi;
% 			else
% 				mask = true(size(xi)); 
% 			end
% 			xi = xi(mask); d = d(mask);
% 			fprintf('%s, %s bw: %.04f\n', ff, mm, bw);
% 			plot(ax.(ff), xi, d, 'color', col, ...
% 				'linestyle', ls, 'Displayname', mm, 'linewidth', 3);
% % 			[N, edges] = histcounts(summ.(ff).(mm), 'normalization', 'pdf');
% % 			plot(ax.(ff), edges(2:end) - diff(edges)/2, N)
% 			hold(ax.(ff), 'on');
% 		catch ME
% 			if ~strcmpi(ME.identifier, 'stats:mvksdensity:NanX')
% 				rethrow(ME)
% 			end
% 		end
% 	end
% 	
% 	hold(ax.(ff), 'off')
% 	legend(ax.(ff))
% 	title(ax.(ff), ff)
% 	ylabel('PDE')
% 	switch ff
% 		case "N_detections_early"
% 			xlabel(ax.(ff), 'Count');
% 			ylabel(ax.(ff), 'Early Detections');
% 		case {"Phi_mean_early", "Phi_mean_late" "Phi_std_early" "Phi_std_late"}
% 			xlabel(ax.(ff), 'Direction');
% 			ylabel(ax.(ff), 'Early Detections');
% 		case "First_detection"
% 			xlabel(ax.(ff), 'Time (s)');
% 	end
% end
% 
% end
% 
% function export_efg_(ax)
% % s = load('figs/bos_v_ny/fig1/fig1_efg_style');
% H = figure('units', 'inches', 'position', [0 0 2.15 1.5]);
% for ff = string(fieldnames(ax)')
% 	
% 	aa = ax.(ff);
% 	h = aa.Parent;
% 	axis(aa, 'tight');
% 	fname = ['figs/bos_v_ny/fig1/f1_' char(ff)];
% 	savefig(h, fname); 
% 	lgd = findobj(h, 'type', 'legend');
% 	lgd.Visible = 'off';
% 	a2 = copyobj(aa, H);
% 	set(a2.Title, 'visible', 'off');
% 	set(a2, 'fontname', 'nimbus sans', 'fontsize', 9, 'fontweight', 'bold', ...
% 		'outerposition', [0 0 1 1], 'linewidth', 2)
% 	print(H, fname, '-dpng');
% 	delete(a2);
% % 	hgexport(h, fname, s);
% end
% end
% 
% function runsplit2fit()
% % convert files from run_split	
% sims = 101:200;
% % sims = 165:198;
% directory_tag = '_d10';
% metric_name = 'D10';
% for ss = sims
% 	files = dir(['FS_Seizure' num2str(ss) directory_tag '/*mat']);
% 	clear temp
% 	varname = string(who('-file', [files(1).folder filesep files(1).name]));
% 	for ii = length(files):-1:1
% 		contents = load([files(ii).folder filesep files(ii).name]);
% 		temp(ii) = contents.(varname);
% 	end
% 	clear out
% 	out.(metric_name) = WaveProp.resize_obj(temp);
% 	save(['FS_Seizure' num2str(ss) '_fits.mat'], '-struct', 'out', '-append');
% end
% end
% 
% function style = set_style_()
% cmap = lines(7);
% style.M.color = cmap(1, :);
% style.E.color = cmap(2, :);
% style.D10.color = cmap(3, :);
% style.D1.color = cmap(5, :);
% style.D1xwh.color = cmap(4, :);
% 
% end
% 
%%
% function save_BCDH_(F, sub)
% 
% assert(ismember(sub, 'BCDH'));
% h = F.(sub);
% aa = h.Children;
% ylabel(aa, 'Direction');
% xlabel(aa, 'Time (s)');
% % ttl_str = struct( ...
% %     'B', F.MetricNames.M, ...
% %     'C', F.MetricNames.E, ...
% %     'D', F.MetricNames.D10, ...
% %     'H', F.MetricNames.D1xwh);
% 
% h_new = figure('position', [0 0 2.5 1.55]);
% aa_new = copyobj(aa, h_new);
% title(aa_new, F.MetricNames.(h.Tag));
% aa_new.OuterPosition = [0 0 1 1];
% 
% fname = [F.prefix 'bcde_' h.Tag];
% print(h_new, fname, '-dpng');
% savefig(h, fname);
% close(h_new);
% 
% end
