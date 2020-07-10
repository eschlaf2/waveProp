%
% Computes averaged phase of wavelet time-frequency decomposition from a
% persyst file. locked on event codes form a .evt file
% 
% Julien Besle 01/XX/2008

function wavelet_phase = wavelet_phase_persyst(event_file, event_codes, low_frequency, high_frequency, start_latency, stop_latency)

wavelet_phase = [];

if nargin ~= 4 && nargin ~= 6
    fprintf (1, 'Usage: wavelet_phase_persyst ( event_file, [vector of event_codess], low_frequency, high_frequency, start_latency (in ms), stop latency(in ms))\n');
    return;
end

% if nargin == 4
%     start_baseline = start_latency;
%     stop_baseline = 0;
% end

omega = 6; %size of the wavelet

%reading persyst file headers
load('-mat',event_file); %contains triggers and data_file
if exist('persyst_triggers','var')
    triggers = persyst_triggers;
    data_file = persyst_lay_file;
end
for i_session = 1:size(triggers,3)
    fileinfo(i_session) = read_persyst (data_file(i_session,:));
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

number_of_samples = (stop_latency - start_latency)*samp_rate/1000;

if low_frequency < samp_rate/number_of_samples
    fprintf(1,'low_frequency must be at least %.0f / %d = %.2f Hz\n', samp_rate, number_of_samples, samp_rate/number_of_samples);
    return
end

if high_frequency > samp_rate/2
    fprintf(1,'high_frequency should not be larger than %.0f/2 = %.1f Hz\n',samp_rate, samp_rate/2);
    return
end


wavelet_phase.samp_rate = samp_rate;
wavelet_phase.zero = (-1) * round (start_latency * samp_rate/1000) + 1;
wavelet_phase.event_codes = event_codes;
wavelet_phase.event_names = event_names(event_codes);
wavelet_phase.Xunit = 'ms';
wavelet_phase.Yunit = 'Hz';
%wavelet_phase.baseline_start = (start_baseline-start_latency)*wavelet_phase.samp_rate/1000+1;
%wavelet_phase.baseline_stop = (stop_baseline-start_latency)*wavelet_phase.samp_rate/1000+1;
wavelet_phase.name = [event_file ' : Phase of wavelet decomposition'];
wavelet_phase.chanlist = fileinfo(1).chanlist;

for i_event = 1:length(event_codes) %for each type of event
    waitbar_string = ['Computing amplitude of wavelet for event ',num2str(event_codes(i_event))];
    h_waitbar = waitbar(0,waitbar_string);
    sum = zeros(fileinfo(1).nchannels, number_of_samples);
    total_trials(i_event) = 0;
%    figure;
%    hold on;
    for i_session = 1:size(triggers,3) %if the trigger file refers to different data files
        trials_times = triggers(1,find(triggers(2,:,i_session)==event_codes(i_event) & triggers(3,:,i_session)>0),i_session); 
        %all the trials with this event number and that have not been excluded
        for i_trials = 1:length(trials_times)
            waitbar(i_trials/length(trials_times),h_waitbar);
            %read the trial form start_latency to stop_latency
            eeg = read_persyst_eeg(fileinfo(i_session), trials_times(i_trials)/samp_rate+start_latency/1000, (stop_latency - start_latency)/1000);
            if size(eeg.eeg_data,2) == number_of_samples %we want only complete trials
                for i_channel = 1:size(eeg.eeg_data,1)
                    [wave,period,scale,cone_of_influence] = basewave4(double(eeg.eeg_data(i_channel,:))',samp_rate,low_frequency,high_frequency,omega,0);
                    wave_phase(i_channel,:,:) = angle(wave); 
                end
                if total_trials(i_event) == 0
                    sum = zeros(size(wave_phase));
                end

                sum = sum + wave;%compute the sum
                total_trials(i_event) = total_trials(i_event) + 1;% the nmber of trials in the sum
            end
            clear eeg.eeg_data
%            plot(wavelet_phase.data(71,:,i_event),'k');
        end
%        fprintf(1,'\n\n');
    end
    
    %divide by the number of trials 
    wavelet_phase.data(:,:,:,i_event) = sum / total_trials(i_event);
    %hold off
    close(h_waitbar);
end
wavelet_phase.number_of_trials = total_trials;
wavelet_phase.cone_of_influence = 1./cone_of_influence;
wavelet_phase.frequencies = 1./period;
