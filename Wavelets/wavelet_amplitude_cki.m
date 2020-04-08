%
% Averaged amplitude of wavelet time-frequency decomposition
% 

function wavelet_amplitude = wavelet_amplitude_cki(event_file, event_codes, low_frequency, high_frequency, start_latency, stop_latency, channel_name)

wavelet_amplitude = [];

if nargin ~= 6 & nargin~=7
    fprintf (1, 'Usage: wavelet_amplitude_cki ( event_file, [vector of event_codess], low_frequency, high_frequency, start_latency (in ms), stop latency(in ms), [channel_name])\n');
    return;
end

% if nargin == 4
%     start_baseline = start_latency;
%     stop_baseline = 0;
% end


omega = 6; %size of the wavelet

%reading persyst file headers
load('-mat',event_file); %contains triggers and data_file

for i_session = 1:size(triggers,3)
    fileinfo(i_session) = read_ns5_header (data_file(i_session,:));
    if i_session == 1
	samp_rate = fileinfo(i_session).samp_rate;
    else
        if fileinfo(i_session).nchannels ~= fileinfo(i_session - 1).nchannels
            fprintf(1,'different sessions must have the same number of channels\n');
            return;
        end
        if fileinfo(i_session).samp_rate ~= fileinfo(1).samp_rate
            fprintf(1,'different sessions must have the same sample rate\n');
            return;
        end
    end
end

% if nargin ==7
%     channel_num = find_index(channel_name,fileinfo(1).chanlist);
%     if ~channel_num
%         fprintf(1,'channel not found\n');
%         return;
%     end
% end
number_of_samples = (stop_latency - start_latency)*samp_rate/1000;

if low_frequency < samp_rate/number_of_samples
    fprintf(1,'low_frequency must be at least %.0f / %d = %.2f Hz\n', samp_rate, number_of_samples, samp_rate/number_of_samples);
    return
end

if high_frequency > samp_rate/2
    fprintf(1,'high_frequency should not be larger than %.0f/2 = %.1f Hz\n',samp_rate, samp_rate/2);
    return
end

% downsample to 4 times high frequency
%downsamp_rate = min (high_frequency * 4, samp_rate);
%decim_rate = round(samp_rate / downsamp_rate);
decim_rate = 1;
downsamp_rate = samp_rate/decim_rate;
reduced_samples = floor(number_of_samples/decim_rate);

wavelet_amplitude.samp_rate = downsamp_rate;
wavelet_amplitude.zero = (-1) * round (start_latency * downsamp_rate / 1000) + 1;
wavelet_amplitude.event_codes = event_codes;
wavelet_amplitude.event_names = event_names(event_codes);
wavelet_amplitude.Xunit = 'ms';
wavelet_amplitude.Yunit = 'Hz';
%wavelet_amplitude.baseline_start = (start_baseline-start_latency)*wavelet_amplitude.samp_rate/1000+1;
%wavelet_amplitude.baseline_stop = (stop_baseline-start_latency)*wavelet_amplitude.samp_rate/1000+1;
wavelet_amplitude.name = [event_file ' : Amplitude of wavelet decomposition'];
wavelet_amplitude.chanlist = fileinfo(1).chanlist;

for i_event = 1:length(event_codes) %for each type of event
    waitbar_string = ['Computing amplitude of wavelet for event ',num2str(event_codes(i_event))];
    h_waitbar = waitbar(0,waitbar_string);
    set(h_waitbar,'Name',data_file(1,:));
    sum = zeros(fileinfo(1).nchannels, reduced_samples);
    reduced_data = zeros(fileinfo(1).nchannels, reduced_samples);
    total_trials(i_event) = 0;
%    figure;
%    hold on;
    for i_session = 1:size(triggers,3) %if the trigger file refers to different data files
        trials_times = triggers(1,find(triggers(2,:,i_session)==event_codes(i_event) & triggers(3,:,i_session)>0),i_session); 
        %all the trials with this event number and that have not been excluded
        for i_trials = 1:length(trials_times)
            waitbar(i_trials/length(trials_times),h_waitbar);
            %read the trial form start_latency to stop_latency
            if nargin == 6
                eeg = read_cki(data_file(i_session,:), trials_times(i_trials)/samp_rate+start_latency/1000, (stop_latency - start_latency)/1000);
            elseif nargin == 7
                eeg = read_cki(data_file(i_session,:), trials_times(i_trials)/samp_rate+start_latency/1000, (stop_latency - start_latency)/1000,channel_name);
                eeg.eeg_data = eeg.eeg_data';
            end
            if size(eeg.eeg_data,2) == number_of_samples %we want only complete trials
                for i_channel = 1:size(eeg.eeg_data,1)
                    try
                    reduced_data(i_channel,:) = decimate (eeg.eeg_data(i_channel,:), decim_rate);
                    catch
                        size(decimate (eeg.eeg_data(i_channel,:), decim_rate))
                        size(reduced_data(i_channel,:))
                        pause;
                    end
                    if high_frequency <= 150
                        [wave,period,scale,cone_of_influence] = basewave4(reduced_data(i_channel,:)',downsamp_rate,low_frequency,high_frequency,omega,0);
                    else
                        [wave,period,scale,cone_of_influence] = HFbasewave4(reduced_data(i_channel,:)',downsamp_rate,low_frequency,high_frequency,omega,0);
                    end
                    %Peter Lakatos's correction
                    wave =(((abs(wave)/24).^2)*1000)/downsamp_rate; 
                    for i_freq=1:size(wave,1)
                        corrected_wave(i_channel,i_freq,:)=sqrt(wave(i_freq,:)./period(i_freq));
                    end
                end
                if total_trials(i_event) == 0
                    sum = zeros(size(corrected_wave));
                end

                sum = sum + corrected_wave;%compute the sum
                total_trials(i_event) = total_trials(i_event) + 1;% the nmber of trials in the sum
            end
            clear eeg.eeg_data
%            plot(wavelet_amplitude.data(71,:,i_event),'k');
        end
%        fprintf(1,'\n\n');
    end
    
    %divide by the number of trials and baseline correct
    %wavelet_amplitude.data(:,:,i_event) = (sum - repmat(mean(sum(:,wavelet_amplitude.baseline_start:wavelet_amplitude.baseline_stop),2),1,number_of_samples)) / total_trials(i_event);
    wavelet_amplitude.data(:,:,:,i_event) = sum / total_trials(i_event);
    %hold off
    close(h_waitbar);
end
wavelet_amplitude.number_of_trials = total_trials;
wavelet_amplitude.cone_of_influence = 1./cone_of_influence;
wavelet_amplitude.frequencies = 1./period;
