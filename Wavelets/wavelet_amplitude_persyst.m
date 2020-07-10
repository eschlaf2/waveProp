%
% Averaged amplitude of wavelet time-frequency decomposition
% 

function [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = wavelet_amplitude_persyst(event_file, event_codes, low_frequency, high_frequency, start_latency, stop_latency, montage, reference_map, analyzed_map)

wavelet_amp = [];
wavelet_resultant = [];
wavelet_plv = [];
p_rayleigh = [];

if  ~ismember(nargin,[6 7 8 9])
    fprintf (1, 'Usage: [wavelet_amp, wavelet_plv, p_rayleigh, wavelet_resultant_phase] = ...\n');
    fprintf (1, 'wavelet_amplitude_persyst ( event_file, [vector of event_codess], low_frequency, high_frequency, start_latency (in ms), stop latency(in ms),...\n');
    fprintf (1, '[<2D derivative|average>], [<reference_map|reference_channels>], [<analyzed_map|analyzed_channels>])\n');
    return;
end

omega = 6; %size of the wavelet
out_samp_rate=500;

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

fprintf(1,'Input data at %d Hz\n',samp_rate);
fprintf(1,'Output wavelet at %d Hz\n',out_samp_rate);
if out_samp_rate> samp_rate
    fprintf(1,'Decrease out_samp_rate');
    return;
end
number_of_samples = (stop_latency - start_latency)*samp_rate/1000;

if low_frequency <= samp_rate/number_of_samples
    fprintf(1,'low_frequency must be more than %.0f / %d = %.2f Hz\n', samp_rate, number_of_samples, samp_rate/number_of_samples);
    return
end

if high_frequency > samp_rate/2
    fprintf(1,'high_frequency should not be larger than %.0f/2 = %.1f Hz\n',samp_rate, samp_rate/2);
    return
end


wavelet_amp.samp_rate = samp_rate/round(samp_rate/out_samp_rate);
wavelet_amp.zero = (-1) * round (start_latency * samp_rate/round(samp_rate/out_samp_rate) / 1000) + 1;
wavelet_amp.event_codes = event_codes;
wavelet_amp.event_names = event_names(event_codes);
wavelet_amp.Xunit = 'ms';
wavelet_amp.Yunit = 'Hz';
%wavelet_amp.baseline_start = (start_baseline-start_latency)*wavelet_amp.samp_rate/1000+1;
%wavelet_amp.baseline_stop = (stop_baseline-start_latency)*wavelet_amp.samp_rate/1000+1;
wavelet_amp.name = '';

%-------------------- List of analyzed channels
if ~exist('analyzed_map','var') || isempty(analyzed_map)
    analyzed_channels = [1:length(ns5_info.chanlist)];
elseif length(analyzed_map)<5 || ~strcmp(analyzed_map(end-3:end),'.map')
    analyzed_channels = find_index(analyzed_map,ns5_info.chanlist);
    if ~analyzed_channels
        fprintf(1,'channel not found\n');
        return;
    else
    end
else
    analyzed_channels = good_channels(analyzed_map);
end
wavelet_amp.chanlist = ns5_info.chanlist(analyzed_channels);

%------------- Referencing
if exist('montage','var') && ~isempty(montage)
    switch(montage)
        case '2D derivative'
            if  (~exist('reference_map','var') || isempty(reference_map))|| length(reference_map)<5 || ~strcmp(reference_map(end-3:end),'.map')
                [FILENAME, PATHNAME] = uigetfile([pwd '/*.map'], 'input *.map file');
                map_filename=[PATHNAME FILENAME];
            end
            [channels,coeffs] = map_channel_combination(reference_map);
            wavelet_amp.name = [wavelet_amp.name ' : 2D derivative (' reference_map ')'];
            
        case 'average'
            if ~exist('reference_map','var') || isempty(reference_map)
                reference_channels = setdiff(ns5_info.chanlist,{'Stim' 'stim' 'Event'});
                wavelet_amp.name = [wavelet_amp.name ' : average reference'];
            elseif iscell(reference_map) || length(reference_map)<5 || ~strcmp(reference_map(end-3:end),'.map')
                reference_channels = find_index(reference_map,ns5_info.chanlist);
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
            channels = [[1:length(reference_channels)]' repmat([1:length(reference_channels)],length(reference_channels),1)];
            coeffs = [ones(length(reference_channels),1) -1/length(reference_channels)*ones(length(reference_channels),length(reference_channels))];
    end
else
    montage = '';
    fprintf(1,'Monopolar montage\n');
end

%------------------------ Other output structures
wavelet_plv = wavelet_amp;
wavelet_plv.Yunit = '';
wavelet_resultant_phase = wavelet_amp;
p_rayleigh = wavelet_plv;
wavelet_amp.name = [event_file ' : Amplitude of wavelet decomposition' wavelet_amp.name];
wavelet_plv.name = [event_file ' : Phase locking factor (Wavelet)' wavelet_plv.name];
wavelet_resultant_phase.name = [event_file ' : Resultant Phase across trials(Wavelet)' wavelet_resultant_phase.name];
wavelet_resultant_phase.Yunit = 'rad';
p_rayleigh.name = [event_file ' : p value of the rayleigh test (Wavelet)' p_rayleigh.name]; 


for i_event = 1:length(event_codes) %for each type of event
    waitbar_string = ['Computing amplitude of wavelet for event ',num2str(event_codes(i_event))];
    h_waitbar = waitbar(0,waitbar_string);
    set(h_waitbar,'Name',data_file(1,:));
    sum = zeros(fileinfo(1).nchannels, number_of_samples);
    total_trials(i_event) = 0;
%    figure;
%    hold on;
    for i_session = 1:size(triggers,3) %if the trigger file refers to different data files
        trials_times = triggers(1,find(triggers(2,:,i_session)==event_codes(i_event) & triggers(3,:,i_session)>0),i_session); 
        %all the trials with this event number and that have not been excluded
        for i_trials = 1:length(trials_times)
            waitbar(i_trials/length(trials_times),h_waitbar);
            %read the trial from start_latency to stop_latency
            if length(used_channels)>1
                eeg = read_persyst_eeg(fileinfo(i_session), trials_times(i_trials)/samp_rate+start_latency/1000, (stop_latency - start_latency)/1000);
                eeg.eeg_data = eeg.eeg_data(used_channels,:);
            else 
                eeg = read_persyst_eeg(fileinfo(i_session), trials_times(i_trials)/samp_rate+start_latency/1000, (stop_latency - start_latency)/1000,wavelet_amp.chanlist{1});
                eeg.eeg_data = eeg.eeg_data';
            end
            if size(eeg.eeg_data,2) == number_of_samples %we want only complete trials
                if ~isempty(montage)
                    data = eeg.eeg_data;
                    for i_channel = 1:size(channels,1)
                        data(channels(i_channel,1),:) = sum(double( eeg.eeg_data(nonzeros(channels(i_channel,:)),:) ...
                            .* repmat(nonzeros(coeffs(i_channel,:)),1,size(eeg.eeg_data,2)) ),1);
                    end
                    eeg.eeg_data = data;
                end
                for i_channel = 1:length(analyzed_channels)
                    [wave,period,scale,cone_of_influence] = basewave4(double(eeg.eeg_data(i_channel,:))',samp_rate,low_frequency,high_frequency,omega,0);
                    wave_phase(i_channel,:,:) = exp(angle(wave)* i); 

                    % % Peter's correction 
                    % flattens the spectrum  and gives a constant outcome when  varying the sample_rate (for a given number of samples), 
                    % corrected_wave(i_channel,:,:) = (abs(wave)/24)./ repmat(sqrt(period)',1,size(wave,2)) * sqrt(1000/samp_rate);

                    %my own correction so that 
                    %   1- the wavelet spectrum has the same scale as the the spectrum before multiplying by f
                    %   2- multiply by sqrt(f) so that the spectrum of a 1/f signal(pink noise) is flat
                    corrected_wave(i_channel,:,:) = abs(wave)*sqrt(5)*sqrt(size(wave,2))./repmat(sqrt(period)',1,size(wave,2));
                    % works fine for omega = 6 or 7
                    
                end
                if total_trials(i_event) == 0
                    sum = zeros(size(corrected_wave));
                    resultant = zeros(size(corrected_wave));
                end

                sum = sum + corrected_wave;%compute the sum
                resultant = resultant + wave_phase;
                total_trials(i_event) = total_trials(i_event) + 1;% the number of trials in the sum
            end
            clear eeg.eeg_data
%            plot(wavelet_amp.data(71,:,i_event),'k');
        end
%        fprintf(1,'\n\n');
    end
    
    %divide by the number of trials and baseline correct
    %wavelet_amp.data(:,:,i_event) = (sum - repmat(mean(sum(:,wavelet_amp.baseline_start:wavelet_amp.baseline_stop),2),1,number_of_samples)) / total_trials(i_event);
    waitbar_string = ['Downsampling data for event ',num2str(event_codes(i_event))];
    waitbar(0, h_waitbar, waitbar_string);
    for i_channel=1:length(analyzed_channels)
        wavelet_amp.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(sum(i_channel,:,:) / total_trials(i_event))',round(samp_rate/out_samp_rate))';
        wavelet_resultant.data(i_channel,:,:,i_event) = averaged_downsample(squeeze(abs(resultant(i_channel,:,:)) / total_trials(i_event))',round(samp_rate/out_samp_rate))';
        waitbar(i_channel/size(eeg.eeg_data,1),h_waitbar);
    end
    %hold off
    close(h_waitbar);
    n = total_trials(i_event);
    Z = n*wavelet_resultant.data(:,:,:,i_event).^2;
    p_rayleigh.data(:,:,:,i_event) = exp(-Z) .* (1 + (2*Z - Z.^2) / (4*n) - (24*Z - 132*Z.^2 + 76*Z.^3 - 9*Z.^4) / (288*n^2));

end
wavelet_amp.number_of_trials = total_trials;
wavelet_resultant.number_of_trials = total_trials;
p_rayleigh.number_of_trials = total_trials;
wavelet_amp.cone_of_influence = 1./averaged_downsample(cone_of_influence,round(samp_rate/out_samp_rate));
wavelet_resultant.cone_of_influence = 1./averaged_downsample(cone_of_influence,round(samp_rate/out_samp_rate));
p_rayleigh.cone_of_influence = 1./averaged_downsample(cone_of_influence,round(samp_rate/out_samp_rate));
wavelet_amp.frequencies = 1./period;
wavelet_resultant.frequencies = 1./period;
p_rayleigh.frequencies = 1./period;
