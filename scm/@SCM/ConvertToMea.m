function ConvertToMea(scm)
    if ~scm.save, return, end
	if ~ismember('label', fieldnames(scm)) || isempty(scm.label), scm.label = 'SCM'; end
    
    TESTING = false;  % This assumes some extra fields are saved in NP
    
	files = dir(sprintf('%s_%d_*mat', scm.basename, scm.sim_num));
	addpath(files(1).folder);
	load(files(1).name, 'last');
	cmap = 1-gray;
	im = round(rescale(last.Ve) * (length(cmap) - 1)) + 1;
	movQ(numel(files) - 1) = im2frame(im, cmap);
	movV(numel(files) - 1) = im2frame(im, cmap);
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
        movK(ii) = im2frame(round(rescale(last.K, 0, 1, 'InputMin', 0, 'InputMax', 1.5) * (length(cmap) - 1)) + 1, cmap);
        tt{ii} = time;
		qe{ii} = NP.Qe;
		ve{ii} = NP.Ve;
        if TESTING
            k{ii} = NP.K;
            qi{ii} = NP.Qi;

            dii{ii} = NP.Dii;
            vi{ii} = NP.Vi;
        end
        ii = ii + 1;
	end
	
	
    qe_mat = cat(1, qe{:});
    if TESTING
        ve_mat = cat(1, ve{:});
        qi_mat = cat(1, qi{:});  % Looking at depolarization block (dbstop here and then <assignin('base', 'vi_mat', vi_mat);>
        k_mat = cat(1, k{:});
        dii_mat = cat(1, dii{:});
        vi_mat = cat(1, vi{:});

        % Just for testing... remove these later
        assignin('base', 've_mat', ve_mat);
        assignin('base', 'qi_mat', qi_mat);
        assignin('base', 'vi_mat', vi_mat);
        assignin('base', 'k_mat', k_mat);
        assignin('base', 'qe_mat', qe_mat);
        assignin('base', 'dii_mat', dii_mat);
    end
	time = cat(1, tt{:});
	sample_rate = min(round(1/mean(diff(time))/1e3)*1e3, scm.subsample);
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
		'Padding', scm.padding, ...
		'Name', [scm.label ' Seizure ' num2str(scm.sim_num)], ...
		'Time', time, ... 
		'Path', scm.mea_path ...	 
		);
	
%   mea = add_noise_(mea, 2);  % Add 3D brownian noise with snr=2; the spectra after this transformation looked similar to recorded seizures - could also use (much) higher snr pink noise
% 	qe_mat = rescale(single(qe_mat), 0, 25);  % range is based on experimentation
% 	mea.firing_rate = reshape(qe_mat, size(mea.Data));
	mea.event_inds = rate2events_(mea);
	mea.event_mat_size = size(mea.Data);
	mea.params = init_mea_params();
	fprintf('Saving %s ... ', mea.Path);
	save(mea.Path, '-struct', 'mea');
	m = matfile(sprintf('%s_%d_info', scm.basename, scm.sim_num), 'Writable', true);
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
