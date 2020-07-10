%
% Averaged amplitude of wavelet time-frequency decomposition
% 

function [p_rand, p_rand_max_sum_stat, p_rand_max_contiguous, p_rand_max_stat, wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase, significant_max_sum_clusters] = ...
    wavelet_permutation_persyst(event_filename, event_codes, low_frequency, high_frequency, start_latency, stop_latency, ...
    mgrid_filename, n_rand, start_analysis, stop_analysis, design, tail, analyzed_map, baseline_start, baseline_stop, mask_tf, event_mask_code)

p_rand = [];
wavelet_amp = [];
wavelet_resultant_phase = [];
wavelet_plv = [];
p_rayleigh = [];

if  ~ismember(nargin,[7 8 10 11 12 13 15 17])
    fprintf (1, 'Usage: [p_rand, p_rand_max_sum_stat, p_rand_max_contiguous, p_rand_max_stat, wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase, significant_max_sum_clusters] = ...\n');
    fprintf (1, '\twavelet_permutation_persyst ( event_filename, [event_codes], low_frequency, high_frequency, start_latency (in ms), stop latency(in ms), ''mgrid_filename'',...\n');
    fprintf (1, '\tn_rand, [start_analysis, stop_analysis], [<''paired''|''independent''], [<''both''|''right''|''left''>], [<''analyzed_map''|}''analyzed_channels''}>],...\n');
    fprintf (1, '\t[baseline_start, baseline_stop], [mask_tf, event_mask_code])\n');
    return;
end

omega = 6; %size of the wavelet
inter_samp_rate = 500;
output_samp_rate=50;
spatial_threshold = 13;  %maximum distance between two contiguous electrodes in mm
single_test_threshold = .05;    %statistical threshold fo the null hypothesis at an individual channel
global_threshold = .05;         %statistical threshold fo the global null hypothesis

if ~exist('design','var') || isempty(design)
    design = 'paired';
end
if ~exist('tail','var') || isempty(tail)
    tail = 'both';
end
if strcmp(design, 'paired') && (~exist('baseline_stop','var') || isempty(baseline_stop))
    baseline_start = start_latency;
    baseline_stop = 0;
end
if exist('mask_tf','var') && ~isempty(mask_tf)
    event_mask_index = find(mask_tf.event_codes==event_mask_code);
    if length(event_mask_index)~=1
        error('Problem with event codes')
    end
end

temp_filename = 'temporary_wavelet_data';
if ~isempty(dir(temp_filename))
    answer = input('Erase temporary_wavelet_data_file (y/n) ?','s');
    if strcmp(answer,'y')
        delete(temp_filename);
    else
        temp_filename = input('Name of temporary data file ?','s');
    end
end

%reading persyst file headers
load('-mat',event_filename); %contains triggers and data_file
if exist('persyst_triggers','var')
    triggers = persyst_triggers;
    data_file = persyst_lay_file;
end
for i_session = 1:size(triggers,3)
    [dump, dump, extension] = fileparts(data_file(i_session,:));
    switch extension
        case '.lay'
            file_info(i_session) = read_persyst (data_file(i_session,:));
    
        case '.eeg'
            file_info(i_session) = read_eeg(data_file(i_session,:),0,1);
    end
    if i_session == 1
        input_samp_rate = file_info(i_session).samp_rate;
    else
        if file_info(i_session).nchannels ~= file_info(i_session - 1).nchannels
            fprintf(1,'different sessions must have the same number of channels\n');
            return;
        end
        if file_info(i_session).samp_rate ~= input_samp_rate
            fprintf(1,'different sessions must have the same sample rate\n');
            return;
        end
    end
end

number_of_samples = (stop_latency - start_latency)*input_samp_rate/1000;
reduc_factor1 = max(1,round(input_samp_rate/inter_samp_rate));
inter_samp_rate = input_samp_rate/reduc_factor1;
reduc_factor2 = max(1,round(inter_samp_rate/output_samp_rate));
output_samp_rate = inter_samp_rate/reduc_factor2;
fprintf(1,'Input data at %d Hz\n',input_samp_rate);
fprintf(1,'processing data at %d Hz\n',inter_samp_rate);
fprintf(1,'Output wavelet at %d Hz\n',output_samp_rate);

if low_frequency <= inter_samp_rate/number_of_samples
    fprintf(1,'low_frequency must be more than %.0f / %d = %.2f Hz\n', inter_samp_rate, number_of_samples, inter_samp_rate/number_of_samples);
    return
end

if high_frequency > inter_samp_rate/2
    fprintf(1,'high_frequency should not be larger than %.0f/2 = %.1f Hz\n',inter_samp_rate, inter_samp_rate/2);
    return
end

if ~exist('stop_analysis','var') || isempty(stop_analysis)
    start_analysis = start_latency;
    stop_analysis = stop_latency;
end
start_analysis = round((start_analysis-start_latency)*output_samp_rate/1000)+1;
stop_analysis = round((stop_analysis-start_latency)*output_samp_rate/1000)+1;

wavelet_amp.samp_rate = output_samp_rate;
wavelet_amp.zero = (-1) * round (start_latency * output_samp_rate / 1000) + 1;
wavelet_amp.event_codes = event_codes;
wavelet_amp.event_names = event_names(event_codes);
wavelet_amp.Xunit = 'ms';
wavelet_amp.Yunit = 'Hz';
wavelet_amp.name = '';

%-------------------- List of analyzed channels
if ~exist('analyzed_map','var') || isempty(analyzed_map)
    analyzed_channels = [1:length(file_info(1).chanlist)];
elseif iscell(analyzed_map) || length(analyzed_map)<5 || ~strcmp(analyzed_map(end-3:end),'.map')
    if ~iscell(analyzed_map)
        analyzed_map = cellstr(analyzed_map);
    end
    analyzed_channels = find_index(analyzed_map,file_info(1).chanlist);
    if length(nonzeros(analyzed_channels))~=length(analyzed_map)
        fprintf(1,'Some channels were not found\n');
        return;
    end
else
    analyzed_channels = good_channels(analyzed_map);
end
wavelet_amp.chanlist = file_info(1).chanlist(analyzed_channels);
n_analyzed_channels = length(analyzed_channels);


%------------------------ Other output structures
wavelet_plv = wavelet_amp;
wavelet_plv.name = [event_filename ' : Phase locking factor (Wavelet)' wavelet_plv.name];
wavelet_resultant_phase = wavelet_amp;
wavelet_resultant_phase.name = [event_filename ' : Resultant Phase across trials(Wavelet)' wavelet_resultant_phase.name];
p_rand = wavelet_plv;
p_rand.name = [event_filename ' : p value of the permutation test (Wavelet)' p_rand.name];
if exist('baseline_stop','var')
    p_rand.baseline_start = round((baseline_start-start_latency)*output_samp_rate/1000)+1;
    p_rand.baseline_stop = round((baseline_stop-start_latency)*output_samp_rate/1000)+1;
end
p_rayleigh = wavelet_plv;
p_rayleigh.name = [event_filename ' : p value of the rayleigh test (Wavelet)' p_rayleigh.name]; 
wavelet_amp.name = [event_filename ' : Amplitude of wavelet decomposition' wavelet_amp.name];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Spatial contiguity%%%%%%%%%%%%%%

%reads the mgrid channel
if ~exist('mgrid_filename','var') || isempty(dir(mgrid_filename))
    [FILENAME, PATHNAME] = uigetfile([pwd '/*.mgrid'], 'input *.mgrid file');
    mgrid_filename = [PATHNAME FILENAME];
end
contiguous = find_contiguous_channels(mgrid_filename,wavelet_amp.chanlist,spatial_threshold);

% for i_channel = 1:length(analyzed_channels)
%     contiguous(i_channel,i_channel)=0;
% end
contiguous = triu(contiguous(1:end-1,:));


%open temporary data file for writing
[temp_file, message] = fopen(temp_filename,'w');
if temp_file <= 0
    error ('Error opening %s for writing: %s', temp_filename, message);
end


h_waitbar = waitbar(0,'','Name',data_file(1,:));
%---------------------------------------- Computation Loop
for i_channel = 1:n_analyzed_channels
    
    for i_event = 1:length(event_codes) %for each type of event
        total_trials(i_event) = 0;
        corrected_wave = [];
        for i_session = 1:size(triggers,3) %if the trigger file refers to different data files
            trials_times = double(triggers(1,find(triggers(2,:,i_session)==event_codes(i_event) & triggers(3,:,i_session)>0),i_session)); 
            %all the trials with this event number and that have not been excluded
            waitbar(0,h_waitbar,['Computing wavelet for channel ',wavelet_amp.chanlist{i_channel}, ', event ' num2str(event_codes(i_event)) ', session ' num2str(i_session)]);
            for i_trials = 1:length(trials_times)
                waitbar(i_trials/length(trials_times),h_waitbar);
                %read the trial from start_latency to stop_latency
                switch extension
                    case '.lay'
                        eeg = read_persyst_eeg(file_info(i_session), trials_times(i_trials)/input_samp_rate+start_latency/1000, (stop_latency - start_latency)/1000, wavelet_amp.chanlist{i_channel});
                        eeg.eeg_data = eeg.eeg_data';
                    case '.eeg'
                        eeg = read_eeg(data_file(i_session,:), trials_times(i_trials)+round(start_latency/1000*input_samp_rate), round((stop_latency - start_latency)/1000*input_samp_rate),analyzed_channels(i_channel));
                end

                if size(eeg.eeg_data,2) == number_of_samples %we want only complete trials
                    total_trials(i_event) = total_trials(i_event) + 1;% the number of trials in the sum_trials
                    %integrate to inter_samp_rate Hz
                    if reduc_factor1>1
                        integrated_data = double(averaged_downsample(eeg.eeg_data',reduc_factor1));
                    else
                        integrated_data = double(eeg.eeg_data)';
                    end
                    [wave,period,scale,cone_of_influence] = basewave4(integrated_data,inter_samp_rate,low_frequency,high_frequency,omega,0);
                    wave_phase = exp(angle(wave)* i); 

                    % % Peter's correction 
                    % flattens the spectrum  and gives a constant outcome when  varying the sample_rate (for a given number of samples), 
                    % corrected_wave(i_channel,:,:) = (abs(wave)/24)./ repmat(sqrt(period)',1,size(wave,2)) * sqrt(1000/inter_samp_rate);

                    %my own correction so that 
                    %   1- the wavelet spectrum has the same scale as the the spectrum before multiplying by f
                    %   2- multiply by sqrt(f) so that the spectrum of a 1/f signal(pink noise) is flat
                    corrected_wave(:,:,total_trials(i_event)) = abs(wave)*sqrt(5)*sqrt(size(wave,2))./repmat(sqrt(period)',1,size(wave,2));
                    % works fine for omega = 6 or 7
%                    corrected_wave(:,:,total_trials(i_event)) = rand(size(corrected_wave,1),size(corrected_wave,2));

                    if total_trials(i_event) == 1
                        sum_trials = zeros(size(corrected_wave,1),size(corrected_wave,2));
                        resultant = sum_trials;
                    end
                    sum_trials = sum_trials + corrected_wave(:,:,total_trials(i_event));%compute the sum_trials
                    resultant = resultant + wave_phase;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Downsampling

        waitbar(0, h_waitbar, ['Downsampling data for channel ',wavelet_amp.chanlist{i_channel}]);
        %output average data (non baseline corrected)
        wavelet_amp.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(sum_trials / total_trials(i_event))',reduc_factor2)';
        wavelet_plv.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(abs(resultant) / total_trials(i_event))',reduc_factor2)';
        wavelet_resultant_phase.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(angle(resultant))',reduc_factor2)';
        %rayleigh test
        n = total_trials(i_event);
        Z = n*wavelet_plv.data(i_channel,:,:,i_event).^2;
        p_rayleigh.data(i_channel,:,:,i_event) = exp(-Z) .* (1 + (2*Z - Z.^2) / (4*n) - (24*Z - 132*Z.^2 + 76*Z.^3 - 9*Z.^4) / (288*n^2));
        
        %data for statistical analysis
        for i_trial=1:size(corrected_wave,3)
            downsampled_data(:,:,i_trial) = averaged_downsample(squeeze(corrected_wave(:,:,i_trial))',reduc_factor2)';
            waitbar(i_trial/size(downsampled_data,3),h_waitbar);
        end
        %baseline correct if paired designed (against baseline) or if requested
        if exist('baseline_stop','var')
            downsampled_data = downsampled_data - repmat(mean(downsampled_data(:,p_rand.baseline_start:p_rand.baseline_stop,:),2),[1 size(downsampled_data,2) 1]);
        end
        
        %reduce to the window of analysis
        data{i_event} = downsampled_data(:,start_analysis:stop_analysis,:);
        clear downsampled_data
%         %%%%%%%% ADD noise for simulations%%%%%%%%%%%%%%%%%%
%         data{i_event} = data{i_event} + randn(size(data{i_event}))*100000;
        
        % saving this to a file for second randomization
        written = fwrite(temp_file,data{i_event},'double');
        if written ~= numel(data{i_event})
            error('Error writing temporary wavelet data file')
        end
    end
    
    
    number_of_frequencies = size(data{1},1);
    number_of_output_samples = size(data{1},2);
    
    if exist('mask_tf','var') && ~isempty(mask_tf)
        if (mask_tf.samp_rate~=output_samp_rate) || (size(mask_tf.data,1) ~= n_analyzed_channels) || (size(mask_tf.data,2) ~= number_of_frequencies) || ~isempty(find((mask_tf.frequencies~=1./period), 1)) 
            error('mask tf not compatible with output')
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% First Randomization %%%%%%%%%%%%%%
    waitbar(0,h_waitbar,['Computing permutations (1/2) for channel ',wavelet_amp.chanlist{i_channel}]);
    
    if i_channel == 1
        switch(design)
            case 'paired'
                for i_event = 1:length(event_codes) %for each type of event
                    permuted_index{i_event}(:,1) = zeros(size(data{i_event},3),1);
                    if size(data{i_event},3)>168
                        for i_rand = 2:n_rand
                            permuted_index{i_event}(:,i_rand) = round(rand(size(data{i_event},3),1));
                        end
                    else
                        %initialize seed
                        mask = create_lfsr_mask(size(data{i_event},3));
                        seed = zeros(size(data{i_event},3),1);
                        while isempty(find(seed, 1))
                            seed = round(rand(size(data{i_event},3),1));
                        end
                        if n_rand > 2^size(data{i_event},3)-1
                            n_rand = 2^size(data{i_event},3)-1;
                            fprintf(1,['Maximum number of permutation = ' str2num(n_rand) '\n']);
                        end
                        permuted_index{i_event}(:,2) = seed;
                        for i_rand = 3:n_rand
                            permuted_index{i_event}(:,i_rand) = next_lfsr(permuted_index{i_event}(:,i_rand-1),mask);
                        end
                    end
                    number_of_tests = length(event_codes);
                    p_rand.number_of_trials = total_trials;
                end
            case 'independent'
                if length(event_codes)==2
                    n_total_trials = total_trials(1) + total_trials(2);
                    number_of_tests = 1;
                    permuted_index{1}(:,1) = 1:n_total_trials;
                    for i_rand = 2:n_rand
                        permuted_index{1}(:,i_rand) = randperm(n_total_trials);
                    end
                    p_rand.number_of_trials = NaN;
                    p_rand.event_codes = 1;
                    p_rand.event_names{1} = [event_names{event_codes(1)} ' - ' event_names{event_codes(2)}];
                else
                    error('Test not implemented for more than two conditions')
                end
              
        end
        p_data = zeros(n_analyzed_channels,number_of_frequencies,number_of_output_samples,number_of_tests);
        threshold = p_data ;
        distribution = zeros(n_rand,number_of_frequencies,number_of_output_samples);
    end
    if strcmp(design,'independent')
        data = {shiftdim([shiftdim(data{1},1) shiftdim(data{2},1)],2)};
    end
    for i_test = 1:number_of_tests
        waitbar(0,h_waitbar,['Computing permutations (1/2) for channel ' wavelet_amp.chanlist{i_channel} ', event ' num2str(event_codes(i_test))]);
        for i_rand = 1:n_rand
            waitbar(i_rand/n_rand,h_waitbar);
            switch(design)
                case 'paired'
                    permuted_data = data{i_test};
                    permuted_data(:,:,permuted_index{i_test}(:,i_rand)==1) = -1*permuted_data(:,:,permuted_index{i_test}(:,i_rand)==1);
                    statistic = squeeze(mean(permuted_data,3));
                case 'independent'
                    if length(event_codes)==2
                        permuted_data = data{i_test}(:,:,permuted_index{i_test}(:,i_rand));
                        statistic = squeeze(mean(permuted_data(:,:,1:total_trials(1)),3)) - squeeze(mean(permuted_data(:,:,total_trials(1)+1:n_total_trials),3));
                    %to modify for 3 events :    
                    %elseif length(event_codes) == 3
                    end
            end

            switch(tail)
                case 'both'
                    statistic = abs(statistic);
                case 'left'
                    statistic = -1*statistic;
            end

            distribution(i_rand,:,:) = statistic;   %register the statistic values  at for this permutation
            p_data(i_channel,:,:,i_test) = p_data(i_channel,:,:,i_test) + (distribution(i_rand,:,:)>distribution(1,:,:)); %count the number of randomization >  the actual value
        end

        sorted_distribution = sort(distribution,1);
        threshold(i_channel,:,:,i_test) = sorted_distribution(round((1-single_test_threshold)*n_rand),:,:); %find the threshold associated with the upper single_test_threshold permutations
        %actual_statistic(i_channel,:,:,i_test) = distribution(1,:,:);
    end
end

fclose(temp_file);

p_data = p_data/(n_rand);

wavelet_amp.number_of_trials = total_trials;
wavelet_plv.number_of_trials = total_trials;
wavelet_resultant_phase.number_of_trials = total_trials;
p_rayleigh.number_of_trials = total_trials;
wavelet_amp.cone_of_influence = 1./averaged_downsample(cone_of_influence,reduc_factor2);
wavelet_plv.cone_of_influence = 1./averaged_downsample(cone_of_influence,reduc_factor2);
wavelet_resultant_phase.cone_of_influence = 1./averaged_downsample(cone_of_influence,reduc_factor2);
p_rayleigh.cone_of_influence = 1./averaged_downsample(cone_of_influence,reduc_factor2);
wavelet_amp.frequencies = 1./period;
wavelet_plv.frequencies = 1./period;
wavelet_resultant_phase.frequencies = 1./period;
p_rayleigh.frequencies = 1./period;

p_rand.frequencies = 1./period;
p_rand.cone_of_influence = wavelet_amp.cone_of_influence(start_analysis:stop_analysis);
if min(p_rand.cone_of_influence) > p_rand.frequencies(1)
    fprintf(1,'Warning: edge effect in analysis window\n');
end
p_rand.zero = wavelet_amp.zero-start_analysis+1;


%other structures output after the second randomization
p_rand_max_stat = p_rand;
p_rand_max_sum_stat = p_rand;
p_rand_max_contiguous = p_rand;
significant_max_sum_clusters = p_rand;

if exist('mask_tf','var') && ~isempty(mask_tf)
    time_shift =  p_rand.zero - mask_tf.zero;
    start_trial_tf = max(1,time_shift+1);
    stop_trial_tf = min(size(p_data,3),size(mask_tf.data,3)+time_shift);
    start_trial_mask = max(1,-time_shift);
    stop_trial_mask = min(size(mask_tf.data,3), size(p_data,3)-time_shift);
    mask_data = mask_tf.data(:,:,start_trial_mask:stop_trial_mask,event_mask_index);
    mask_data = mask_data~=0;
%     figure('name','Mask');
%     hold on;
%     for i_channel = 1:size(mask_data,1)
%         subplot(6,8,i_channel);
%         surface(squeeze(mask_data(i_channel,:,:)+0),'EdgeColor','none');
%     end
%     h_unmasked = figure('name','Unmasked');
%     h_masked = figure('name','Masked');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Second randomization
waitbar(0,h_waitbar,['Computing permutations (2/2) (all channels)']);
p_rand_max_stat.data = ones(n_analyzed_channels,number_of_frequencies,number_of_output_samples,number_of_tests);
actual_statistic = zeros(size(p_data)) ;
statistic = zeros(size(p_data)) ;
max_stat = zeros(n_rand,number_of_tests);
max_sum_stat = zeros(n_rand,number_of_tests);
max_contiguous = zeros(n_rand,number_of_tests);

for i_rand = 1:n_rand
    waitbar(i_rand/n_rand,h_waitbar);

    [temp_file, message] = fopen(temp_filename,'r');
    if temp_file <= 0
        error ('Error opening %s for reading: %s', temp_filename, message);
    end

    for i_channel = 1:n_analyzed_channels
        for i_test = 1:length(event_codes)
            data{i_test} = fread(temp_file,number_of_frequencies*number_of_output_samples*total_trials(i_test),'double');
            data{i_test} = reshape(data{i_test},[number_of_frequencies number_of_output_samples total_trials(i_test)]);
        end
        if strcmp(design,'independent')
            data = {shiftdim([shiftdim(data{1},1) shiftdim(data{2},1)],2)};
        end
    
        for i_test = 1:number_of_tests
            switch(design)
                case 'paired'
                    permuted_data = data{i_test};
                    permuted_data(:,:,permuted_index{i_test}(:,i_rand)==1) = -1*permuted_data(:,:,permuted_index{i_test}(:,i_rand)==1);
                    statistic(i_channel,:,:,i_test) = shiftdim(squeeze(mean(permuted_data,3)),-1);
                case 'independent'
                    if length(event_codes)==2
                        permuted_data = data{i_test}(:,:,permuted_index{i_test}(:,i_rand));
                        statistic(i_channel,:,:,i_test) = shiftdim(squeeze(mean(permuted_data(:,:,1:total_trials(1)),3)) - squeeze(mean(permuted_data(:,:,total_trials(1)+1:n_total_trials),3)),-1);
                    %to modify for 3 events :    
                    %elseif length(event_codes) == 3
                    end
            end

        end
    end
    switch(tail)
        case 'both'
            statistic = abs(statistic);
        case 'left'
            statistic = -1*statistic;
    end

    if i_rand ==1
        actual_statistic = statistic;
    end


    %-----------global statistics : contiguous clusters of significant channels*frequencies*samples    
    thresholded_statistic = statistic>threshold;
    if exist('mask_tf','var') && ~isempty(mask_tf)
        for i_test = 1:number_of_tests
%             figure(h_unmasked);
%             for i_channel = 1:size(thresholded_statistic,1)
%                 subplot(6,8,i_channel);
%                 surface(squeeze(thresholded_statistic(i_channel,:,:,i_test)+0),'EdgeColor','none');
%             end
            masked_data = thresholded_statistic(:,:,start_trial_tf:stop_trial_tf,i_test);
            masked_data(mask_data==0) = 0;
            thresholded_statistic(:,:,:,i_test) = zeros(n_analyzed_channels,number_of_frequencies,number_of_output_samples);
            thresholded_statistic(:,:,start_trial_tf:stop_trial_tf,i_test) = masked_data;
%             figure(h_masked);
%             for i_channel = 1:size(thresholded_statistic,1)
%                 subplot(6,8,i_channel);
%                 surface(squeeze(thresholded_statistic(i_channel,:,:,i_test)+0),'EdgeColor','none');
%             end
            %-----------global statistic = max accross all electrodes*frequencies*samples
            masked_statistic = statistic(:,:,start_trial_tf:stop_trial_tf,i_test);
            max_stat(i_rand,i_test) =  squeeze(max(max(max(masked_statistic(mask_data~=0)))));
        end
    else
        %-----------global statistic = max accross all electrodes*frequencies*samples
        max_stat(i_rand,:) =  squeeze(max(max(max(statistic))));
    end


    for i_test = 1:number_of_tests
        temp_statistic = statistic(:,:,:,i_test);
        try
         [clustered_statistic,cluster_size] = find_3D_clusters(thresholded_statistic(:,:,:,i_test), contiguous);
        catch
            fprintf('Problem with find_3D_clusters.m at randomization #%d\n',i_rand);
            break;
        end
        max_contiguous(i_rand,i_test) = max(cluster_size);
        
        for i_cluster = 1:length(cluster_size)
            cluster_sum = sum(temp_statistic(clustered_statistic==i_cluster));
            if cluster_sum > max_sum_stat(i_rand,i_test)
                max_sum_stat(i_rand,i_test) = cluster_sum;
            end
            if i_rand ==1
                actual_cluster_sum{i_test}(i_cluster) = cluster_sum;
            end
        end
        
        if i_rand ==1
            actual_clustered_statistic(:,:,:,i_test) = clustered_statistic;
            actual_cluster_size{i_test} = cluster_size;
        end
    end
    
    fclose(temp_file);
end

%-----------global statistic = max accross all electrodes*frequencies*samples
sorted_max_stat = sort(max_stat,1); %order the global statistic max_stat accross permutations
threshold_max_stat = repmat(shiftdim(sorted_max_stat(ceil(n_rand*(1-global_threshold)),:),-2),[n_analyzed_channels number_of_frequencies number_of_output_samples 1]);
p_rand_max_stat.data(actual_statistic>threshold_max_stat) = p_data(actual_statistic>threshold_max_stat);
p_rand_max_stat.name = [event_filename ' : p value of the permutation test (Wavelet) corrected with max(stat) (p = ' num2str(global_threshold) ')'];

%-----------global statistic : max(sum(statistic)) of significant clusters  across channels*frequencies*samples    
%-----------global statistic : max of number of contiguous significant channels*frequencies*samples    
sorted_max_contiguous = sort(max_contiguous,1); %order the global statistic max_contiguous accross permutations
sorted_max_sum_stat = sort(max_sum_stat,1); %order the global statistic max_sum_stat accross permutations
threshold_max_sum_stat = sorted_max_sum_stat(ceil(n_rand*(1-global_threshold)),:);
threshold_max_contiguous = sorted_max_contiguous(ceil(n_rand*(1-global_threshold)),:);

for i_test = 1:number_of_tests
    significant_clusters = 0;
    temp_statistic = actual_clustered_statistic(:,:,:,i_test);
    temp_data = p_data(:,:,:,i_test);
    p_data_max_sum_stat = ones(size(temp_statistic));  %this has size n_analyzed_channels*number_of_frequencies*number_of_output_samples
    p_data_clusters = zeros(size(temp_statistic));  %this has size n_analyzed_channels*number_of_frequencies*number_of_output_samples
    p_data_max_contiguous = ones(size(temp_statistic));
    for i_cluster = 1:length(actual_cluster_size{i_test})
        if actual_cluster_sum{i_test}(i_cluster)>threshold_max_sum_stat(i_test)
            significant_clusters = significant_clusters + 1;
            p_data_clusters(temp_statistic==i_cluster) = significant_clusters;
            p_data_max_sum_stat(temp_statistic==i_cluster) = temp_data(temp_statistic==i_cluster);
        end
        if actual_cluster_size{i_test}(i_cluster)>threshold_max_contiguous(i_test)
            p_data_max_contiguous(temp_statistic==i_cluster) = temp_data(temp_statistic==i_cluster);
        end
    end
    p_rand_max_sum_stat.data(:,:,:,i_test) = p_data_max_sum_stat;
    p_rand_max_contiguous.data(:,:,:,i_test) = p_data_max_contiguous;
    significant_max_sum_clusters.data(:,:,:,i_test) = p_data_clusters;
end

p_rand_max_sum_stat.name = [event_filename ' - p randomization test - corrected with max(sum(stat)) (p = ' num2str(global_threshold) ')'];
p_rand_max_contiguous.name = [event_filename ' - p randomization test - corrected with max(size(cluster)) (p = ' num2str(global_threshold) ')'];
significant_max_sum_clusters.name = [event_filename ' - clusters for max(size(cluster)) correction (p = ' num2str(global_threshold) ')'];
for i_test = 1:number_of_tests
    p_rand_max_stat.event_names{i_test} = [p_rand_max_stat.event_names{i_test} ' [max(stat) = ' num2str(threshold_max_stat(i_test)) ']'];
    p_rand_max_sum_stat.event_names{i_test} = [p_rand_max_sum_stat.event_names{i_test} ' [max(sum(stat)) = ' num2str(threshold_max_sum_stat(i_test)) ']'];
    significant_max_sum_clusters.event_names{i_test} = [significant_max_sum_clusters.event_names{i_test} ' [max(sum(stat)) = ' num2str(threshold_max_sum_stat(i_test)) ']'];
    p_rand_max_contiguous.event_names{i_test} = [p_rand_max_contiguous.event_names{i_test} ' [max(size(cluster)) = ' num2str(threshold_max_contiguous(i_test)) ']'];
end

close(h_waitbar);

delete(temp_filename);

p_rand.data = p_data;   %I put this here, because I don't want duplicates of the data when allocating p_rand_max structures
