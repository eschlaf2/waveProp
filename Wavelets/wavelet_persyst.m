% 
% Averaged amplitude of wavelet time-frequency decomposition
% 

function [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = wavelet_persyst(event_filename, event_codes, low_frequency, high_frequency, start_latency, stop_latency, ...
                                                                           montage, reference_map, analyzed_map)

wavelet_amp = [];
wavelet_resultant_phase = [];
wavelet_plv = [];
p_rayleigh = [];

if  ~ismember(nargin,[6 7 8 9])
    fprintf (1, 'Usage: [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = ...\n');
    fprintf (1, '\twavelet_persyst (''event_filename'', [event_codes], low_frequency, high_frequency, start_latency (in ms), stop latency(in ms),...\n');
    fprintf (1, '\t[<''2D derivative''|''average''>], [<''reference_map''|{''reference_channels''}>], [<''analyzed_map''|{''analyzed_channels''}>])\n');
    return;
end

omega = 6; %size of the wavelet
inter_samp_rate = 5000;
output_samp_rate = 1000;

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

%------------- Referencing
if exist('montage','var') && ~isempty(montage)
    switch(montage)
        case '2D derivative'
            if  (~exist('reference_map','var') || isempty(reference_map))|| length(reference_map)<5 || ~strcmp(reference_map(end-3:end),'.map')
                [FILENAME, PATHNAME] = uigetfile([pwd '/*.map'], 'input *.map file');
                reference_map = [PATHNAME FILENAME];
            end
            [channels,coeffs] = map_channel_combination(reference_map);
            wavelet_amp.name = [wavelet_amp.name ' : 2D derivative (' reference_map ')'];
            
        case 'average'
            if ~exist('reference_map','var') || isempty(reference_map)
                [dump,reference_channels] = setdiff(file_info(1).chanlist,{'Stim' 'stim' 'Event'});
                reference_channels = sort(reference_channels);
                wavelet_amp.name = [wavelet_amp.name ' : average reference (all channels)'];
            elseif iscell(reference_map) || length(reference_map)<5 || ~strcmp(reference_map(end-3:end),'.map')
                reference_channels = find_index(reference_map,file_info(1).chanlist);
                if ~reference_channels
                    error(['Channel not found\n']);
                end
                if iscell(reference_map)
                    ref_string = '';
                    for i_channel = 1:length(reference_map)
                        ref_string = [ref_string ',' reference_map{i_channel}];
                    end
                else
                    ref_string = reference_map;
                end
                wavelet_amp.name = [wavelet_amp.name ' : Referenced to (' ref_string ')'];
            else
                reference_channels = good_channels(reference_map);
                wavelet_amp.name = [wavelet_amp.name ' : average reference (' reference_map ')'];
            end
            channels = [[1:length(file_info(1).chanlist)]' repmat(reference_channels,length(file_info(1).chanlist),1)];
            coeffs = [ones(length(file_info(1).chanlist),1) -1/length(reference_channels)*ones(length(file_info(1).chanlist),length(reference_channels))];

        case 'monopolar'
            montage = '';
            
        otherwise
            error('Montage must be ''2D derivative'', ''average'' or ''monopolar''');
            end
else
    montage = '';
    fprintf(1,'Monopolar montage\n');
end

%------------------------ Other output structures
wavelet_plv = wavelet_amp;
wavelet_plv.name = [event_filename ' : Phase locking factor (Wavelet)' wavelet_plv.name];
wavelet_resultant_phase = wavelet_amp;
wavelet_resultant_phase.name = [event_filename ' : Resultant Phase across trials(Wavelet)' wavelet_resultant_phase.name];
p_rayleigh = wavelet_plv;
p_rayleigh.name = [event_filename ' : p value of the rayleigh test (Wavelet)' p_rayleigh.name]; 
wavelet_amp.name = [event_filename ' : Amplitude of wavelet decomposition' wavelet_amp.name];


%---------------------------------------- Computation Loop
for i_event = 1:length(event_codes) %for each type of event
    waitbar_string = ['Computing wavelet for event ',num2str(event_codes(i_event))];
    h_waitbar = waitbar(0,waitbar_string,'Name',data_file(1,:));
    total_trials(i_event) = 0;
    for i_session = 1:size(triggers,3) %if the trigger file refers to different data files
        trials_times = double(triggers(1,find(triggers(2,:,i_session)==event_codes(i_event) & triggers(3,:,i_session)>0),i_session)); 
        %all the trials with this event number and that have not been excluded
        for i_trials = 1:length(trials_times)
            waitbar(i_trials/length(trials_times),h_waitbar);
            %read the trial from start_latency to stop_latency
            switch extension
                case '.lay'
                    eeg = read_persyst_eeg(file_info(i_session), trials_times(i_trials)/input_samp_rate+start_latency/1000, (stop_latency - start_latency)/1000);

                case '.eeg'
                    eeg = read_eeg(data_file(i_session,:), trials_times(i_trials)+round(start_latency/1000*input_samp_rate), round((stop_latency - start_latency)/1000*input_samp_rate));
            end

            if size(eeg.eeg_data,2) == number_of_samples %we want only complete trials
                %integrate to inter_samp_rate Hz
                if reduc_factor1>1
                    integrated_data = double(averaged_downsample(eeg.eeg_data',reduc_factor1))';
                else
                    integrated_data = double(eeg.eeg_data);
                end
                if ~isempty(montage)
                    data = integrated_data;
                    for i_channel = 1:size(channels,1)
                        data(channels(i_channel,1),:) = sum( integrated_data(nonzeros(channels(i_channel,:)),:) ...
                            .* repmat(nonzeros(coeffs(i_channel,:)),1,size(integrated_data,2)),1);
                    end
                    data = data(analyzed_channels,:);
                else
                    data = integrated_data(analyzed_channels,:);
                end
                
                for i_channel = 1:length(analyzed_channels)
                    [wave,period,scale,cone_of_influence] = basewave4(double(data(i_channel,:))',inter_samp_rate,low_frequency,high_frequency,omega,0);
                    wave_phase(i_channel,:,:) = exp(angle(wave)* i); 

                    % % Peter's correction 
                    % flattens the spectrum  and gives a constant outcome when  varying the sample_rate (for a given number of samples), 
                    % corrected_wave(i_channel,:,:) = (abs(wave)/24)./ repmat(sqrt(period)',1,size(wave,2)) * sqrt(1000/inter_samp_rate);

                    %my own correction so that 
                    %   1- the wavelet spectrum has the same scale as the the spectrum before multiplying by f
                    %   2- multiply by sqrt(f) so that the spectrum of a 1/f signal(pink noise) is flat
                    corrected_wave(i_channel,:,:) = abs(wave)*sqrt(5)*sqrt(size(wave,2))./repmat(sqrt(period)',1,size(wave,2));
                    % works fine for omega = 6 or 7
                    
                end
                if total_trials(i_event) == 0
                    sum_trials = zeros(size(corrected_wave));
                    resultant = sum_trials;
                end

                sum_trials = sum_trials + corrected_wave;%compute the sum_trials
                resultant = resultant + wave_phase;
                total_trials(i_event) = total_trials(i_event) + 1;% the number of trials in the sum_trials
            end
        end
    end
    
    %divide by the number of trials and baseline correct
    %wavelet_amp.data(:,:,i_event) = (sum_trials - repmat(mean(sum_trials(:,wavelet_amp.baseline_start:wavelet_amp.baseline_stop),2),1,number_of_samples)) / total_trials(i_event);
    waitbar_string = ['Downsampling data for event ',num2str(event_codes(i_event))];
    waitbar(0, h_waitbar, waitbar_string);
    for i_channel=1:length(analyzed_channels)
        wavelet_amp.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(sum_trials(i_channel,:,:) / total_trials(i_event))',reduc_factor2)';
        wavelet_plv.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(abs(resultant(i_channel,:,:)) / total_trials(i_event))',reduc_factor2)';
        wavelet_resultant_phase.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(angle(resultant(i_channel,:,:)))',reduc_factor2)';
        waitbar(i_channel/size(data,1),h_waitbar);
    end
    close(h_waitbar);
    n = total_trials(i_event);
    %rayleigh test
    Z = n*wavelet_plv.data(:,:,:,i_event).^2;
    p_rayleigh.data(:,:,:,i_event) = exp(-Z) .* (1 + (2*Z - Z.^2) / (4*n) - (24*Z - 132*Z.^2 + 76*Z.^3 - 9*Z.^4) / (288*n^2));

end

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
