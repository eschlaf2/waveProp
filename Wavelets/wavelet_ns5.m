% computes evoked potentials from a ns5 file, locked to event codes from a .evt file
%
% Julien Besle 04/27/2009

function [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = wavelet_ns5(event_filename, event_codes, low_frequency, high_frequency, start_latency, stop_latency,...
                                                                montage, reference_map, f_design, low_cutting_frq, high_cutting_frq, analyzed_map)

wavelet_amp = [];
wavelet_plv = [];
wavelet_resultant_phase = [];
p_rayleigh = [];

if ~ismember (nargin, [6 7 8 10 11 12])
    fprintf (1, 'Usage: [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = ...\n');
    fprintf (1, '\twavelet_amp_ns5(event_filename, [event_codes], low_frequency, high_frequency,  start_latency (in ms), stop latency(in ms),...\n');
    fprintf (1, '\t[<''bipolar''|''2D derivative''|''average''>], [''reference_map''|{''reference_chanlist''}],...\n');
    fprintf (1, '\t[<''cki noise''|''butter''|''both''>],[low_cutting_frequency],[high_cutting_frequency], [analyzed_map| analyzed_chanlist])\n');
    return;
end


%reading persyst file headers
load('-mat',event_filename); %contains triggers and data_file

omega = 6; %size of the wavelet
inter_samp_rate = 15000;
out_samp_rate = 1000;
file_info = read_ns5_header(data_file);
number_of_samples = (stop_latency - start_latency)*file_info(1).samp_rate/1000;%should add 1 sample ot account for the zero sample but would have to modify the call to read_cki

reduc_factor1 = max(1,round(file_info(1).samp_rate/inter_samp_rate));
inter_samp_rate = file_info(1).samp_rate/reduc_factor1;
reduc_factor2 = max(1,round(inter_samp_rate/out_samp_rate));
out_samp_rate = inter_samp_rate/reduc_factor2;
fprintf(1,'Input data at %d Hz\n',file_info(1).samp_rate);
fprintf(1,'processing data at %d Hz\n',inter_samp_rate);
fprintf(1,'Output wavelet at %d Hz\n',out_samp_rate);

if low_frequency <= inter_samp_rate/number_of_samples
    fprintf(1,'low_frequency must be more than %.0f / %d = %.2f Hz\n', inter_samp_rate, number_of_samples, inter_samp_rate/number_of_samples);
    return
end

if high_frequency > inter_samp_rate/2
    fprintf(1,'high_frequency should not be larger than %.0f/2 = %.1f Hz\n', inter_samp_rate, inter_samp_rate/2);
    return
end


if ~exist('high_cutting_frq','var')
    high_cutting_frq = [];
end
if ~exist('low_cutting_frq','var')
    low_cutting_frq = [];
end
if ~isempty(high_cutting_frq) && inter_samp_rate<high_cutting_frq;
    high_cutting_frq = inter_samp_rate;
end

wavelet_amp.samp_rate = out_samp_rate;
wavelet_amp.zero = (-1) * round (start_latency * out_samp_rate / 1000) + 1;
wavelet_amp.event_codes = event_codes;
wavelet_amp.event_names = event_names(event_codes);
wavelet_amp.Xunit = 'ms';
wavelet_amp.Yunit = 'Hz';
%wavelet_amp.baseline_start = (start_baseline-start_latency)*wavelet_amp.samp_rate/1000+1;
%wavelet_amp.baseline_stop = (stop_baseline-start_latency)*wavelet_amp.samp_rate/1000+1;
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

%----------------------------- Filter design
if exist('f_design','var') && ~isempty(f_design) && ismember(f_design,{'cki noise' 'butter' 'both'})
    if strcmp(f_design,'butter') || strcmp(f_design,'both')
        if ~isempty(low_cutting_frq)
            [b_high,a_high] = filter_design(low_cutting_frq,inter_samp_rate,'high','butter');
            wavelet_amp.name = [wavelet_amp.name ' -  highpass ' num2str(high_cutting_frq) 'Hz'];
        end
        if ~isempty(high_cutting_frq)
            [b_low,a_low] = filter_design(high_cutting_frq,inter_samp_rate,'low','butter');
            wavelet_amp.name = [wavelet_amp.name ' -  lowpass ' num2str(low_cutting_frq) 'Hz'];
        end
    end
    if strcmp(f_design,'cki noise') || strcmp(f_design,'both')
        if isempty(high_cutting_frq)
            cki_cutting_frequencies(2) = inter_samp_rate;
        else
            cki_cutting_frequencies(2) = high_cutting_frq;
        end
        if isempty(low_frequency)
            cki_cutting_frequencies(1) = 0;
        else
            cki_cutting_frequencies(1) = low_cutting_frq;
        end
        wavelet_amp.name = [wavelet_amp.name ' -  cki-noise-filtered'];
    end
else 
    f_design = '';
end

%------------- Referencing
if exist('montage','var') && ~isempty(montage)
    switch(montage)
        case 'bipolar'
            [channels,bipolar_chanlist] = cki_bipolar_pairs(file_info(1).chanlist);
            coeffs = [ones(length(channels),1) -1*ones(length(channels),1)];
            file_info(1).chanlist(channels(:,1)) = bipolar_chanlist;
            wavelet_amp.chanlist = file_info(1).chanlist(analyzed_channels);
            wavelet_amp.name = [wavelet_amp.name ' : bipolar '];

        case '2D derivative'
            if  (~exist('reference_map','var') || isempty(reference_map))|| length(reference_map)<5 || ~strcmp(reference_map(end-3:end),'.map')
                reference_map = '/cki/Patients/maps/cki.map';
                fprintf(1,['Will use ' reference_map ' for 2D derivative computation\n']);
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
            error('Montage must be ''2D derivative'', ''bipolar'',  ''average'' or ''monopolar''');
        
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
            eeg = read_cki(data_file,trials_times(i_trials)/file_info(1).samp_rate+start_latency/1000,(stop_latency - start_latency)/1000);
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
                        data(channels(i_channel,1),:) = sum(double( integrated_data(nonzeros(channels(i_channel,:)),:) ...
                            .* repmat(nonzeros(coeffs(i_channel,:)),1,size(integrated_data,2)) ),1);
                    end
                    integrated_data = data(analyzed_channels,:);
                else
                    integrated_data = integrated_data(analyzed_channels,:);
                end

                if strcmp(f_design,'cki noise') || strcmp(f_design,'both')
                    % CAS modification: add filtering
                    integrated_data = cki_line_noise_filter (integrated_data', inter_samp_rate, cki_cutting_frequencies)';
                end
                
                if strcmp(f_design,'butter') || strcmp(f_design,'both')
                    if ~isempty(high_cutting_frq)
                        integrated_data = filtfilt(b_low,a_low,integrated_data')';
                    end
                    if ~isempty(low_cutting_frq)
                        integrated_data = filtfilt(b_high,a_high,integrated_data')';
                    end
                end

                for i_channel = 1:length(analyzed_channels)
                    [wave,period,scale,cone_of_influence] = basewave4(double(integrated_data(i_channel,:))',inter_samp_rate,low_frequency,high_frequency,omega,0);
                    wave_phase(i_channel,:,:) = exp(angle(wave)* i); 

                    % % Peter's correction 
                    % flattens the spectrum  and gives a constant outcome when  varying the sample_rate (for a given number of samples), 
                    % corrected_wave(i_channel,:,:) = (abs(wave)/24)./ repmat(sqrt(period)',1,size(wave,2)) * sqrt(1000/inter_samp_rate);

                    %my own correction so that 
                    %   1- the wavelet spectrum has the same scale as the fft spectrum before multiplying by f
                    %   2- multiply by f^2 so that the spectrum of a 1/f signal(pink noise) is flat
                    corrected_wave(i_channel,:,:) = abs(wave)*sqrt(5)*sqrt(size(wave,2))./repmat(sqrt(period)',1,size(wave,2));
                    % works fine for omega = 6 or 7

                end
                if total_trials(i_event) == 0
                    sum_trials = zeros(size(corrected_wave));
                    resultant = zeros(size(corrected_wave));
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
        waitbar(i_channel/size(eeg.eeg_data,1),h_waitbar);
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
