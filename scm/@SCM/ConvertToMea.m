function ConvertToMea(PM)
    if ~PM.save, return, end
	if ~isfield(PM, 'label'), PM.label = 'SCM'; end
	files = dir(sprintf('%s_%d_*mat', PM.basename, PM.sim_num));
	addpath(files(1).folder);
	load(files(1).name, 'last');
	cmap = 1-gray;
	im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
	movQ(numel(files) - 1) = im2frame(im, cmap);
	movV(numel(files) - 1) = im2frame(im, cmap);
	[qe, ve, tt] = deal(cell(numel(files) - 1 , 1));
    movK(numel(files) - 1) = im2frame(im, cmap);
	[qe, ve, k, qi, tt, dii, vi] = deal(cell(numel(files) - 1 , 1));

    file_inds = cellfun(@(f) strsplit(f, {'_', '.'}), {files.name}, 'uni', 0);
    file_inds = cellfun(@(f) str2double(strrep(f{end - 1}, 'M', '-')), file_inds);
    [~, file_order] = sort(file_inds, 'ascend');
    
	ii = 1;
    for f = files(file_order)'
		if strfind(f.name, 'info'), continue, end
		load(f.name, 'last', 'NP', 'time');
		disp(f.name)
		disp(ii)
		im = round(rescale(last.Qe) * (length(cmap) - 1)) + 1;
		movQ(ii) = im2frame(im, cmap);
		movV(ii) = im2frame(round(rescale(last.Ve) * (length(cmap) - 1)) + 1, cmap);
        movK(ii) = im2frame(round(rescale(last.K, 0, 1, 'InputMin', 0, 'InputMax', 1) * (length(cmap) - 1)) + 1, cmap);
		qe{ii} = NP.Qe;
		ve{ii} = NP.Ve;
		tt{ii} = time;
        ii = ii + 1;
	end
	
% 	ve_mat = -cat(1, ve{:});
	qe_mat = cat(1, qe{:});
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, PM.subsample);
	dt = 1 / sample_rate;
	nt = size(qe_mat, 1);
	inds = interp1(time, 1:nt, time(1):dt:time(end), 'nearest');
	time =@() time(1):dt:time(end);
% 	ve_mat = ve_mat(inds, :, :);
	qe_mat = qe_mat(inds, :, :);
	
	
	mea = create_mea( ...
		-qe_mat, ... 
        'firing_rate', qe_mat, ...
		'SamplingRate', sample_rate, ... 
		'Padding', PM.padding, ...
		'Name', [PM.label ' Seizure ' num2str(PM.sim_num)], ...
		'Time', time, ... 
		'Path', sprintf('%s/%s/%s_Seizure%d_Neuroport_%d_%d.mat', ...
			pwd, PM.label, PM.label, PM.sim_num, PM.padding) ...	 
		);
	
%   mea = add_noise_(mea, 2);  % Add 3D brownian noise with snr=2; the spectra after this transformation looked similar to recorded seizures - could also use (much) higher snr pink noise
	qe_mat = rescale(single(qe_mat), 0, 25);  % range is based on experimentation
	mea.firingRate = reshape(qe_mat, size(mea.Data));
	mea.event_inds = rate2events_(mea);
	mea.event_mat_size = size(mea.Data);
	mea.params = init_mea_params();
	fprintf('Saving %s ... ', mea.Path);
	save(mea.Path, '-struct', 'mea');
	m = matfile(sprintf('%s_%d_info', PM.basename, PM.sim_num), 'Writable', true);
	m.Qe_movie = movQ;
	m.Ve_movie = movV;
    m.K_movie = movK;
	
	fprintf('Done.\n')
	
	fprintf('Done.\n')
	
end

%% Local functions
function event_inds = rate2events_(mea)
% 	lambda = mea.firingRate;
% 	X = rand(size(lambda));
% 	events = X > exp(-lambda / mea.SamplingRate);
% 	event_inds = find(events);
	
	FR = double(-mea.Data);
	FR = (FR - min(FR)) ./ std(FR);
	[~, event_inds, ~, ~] = findpeaks(FR(:), 'MinPeakHeight', 2);  
end
