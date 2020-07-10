%
% Amplitude of wavelet time-frequency decomposition for one signal
% 

function wavelet = wavelet_amplitude_sig (sig, Fs, low_frequency, high_frequency)

if nargin < 2
    fprintf (1, 'Usage: wavelet_amplitude_persyst (signal, samp_rate, low_frequency, high_frequency\n');
    return;
end
if nargin < 3
    low_frequency = max (1, length(sig)/Fs);
end
if nargin < 4
    high_frequency = Fs/2;
end

omega = 6; %size of the wavelet

number_of_samples = length(sig);
samp_rate = Fs;

if low_frequency < samp_rate/number_of_samples
    warning('low_frequency must be at least %.0f / %d = %.2f Hz\n', samp_rate, number_of_samples, samp_rate/number_of_samples);
    return
end

if high_frequency > samp_rate/2
    warning('high_frequency should not be larger than %.0f/2 = %.1f Hz\n',samp_rate, samp_rate/2);
    return
end


wavelet.samp_rate = samp_rate;
wavelet.zero = 1;
wavelet.Xunit = 'ms';
wavelet.Yunit = 'Hz';
wavelet.name = ['Wavelet decomposition'];

if high_frequency > 200 % looking for high frequency activity
    [wave,period,scale,cone_of_influence] = HFbasewave4(sig',samp_rate,low_frequency,high_frequency,omega,0);
else
    [wave,period,scale,cone_of_influence] = basewave4(sig',samp_rate,low_frequency,high_frequency,omega,0);
end
wavelet.wave = wave;
%wavelet.phase = angle(wave);
%Peter Lakatos's correction
wave =(((abs(wave)/24).^2)*1000)/samp_rate; 
for i_freq=1:size(wave,1)
    % Peter's original correction
    %corrected_wave(i_freq,:)=sqrt(wave(i_freq,:)./period(i_freq));
    % Cathy's correction for 1/f curve - works better for higher
    % frequencies
    %corrected_wave(i_freq,:) = wave(i_freq,:)./sqrt(period(i_freq));
    corrected_wave(i_freq,:) = wave(i_freq,:).*sqrt(i_freq);
end

wavelet.cone_of_influence = 1./cone_of_influence;
wavelet.frequencies = 1./period;
wavelet.corrected_amplitude = corrected_wave;
wavelet.amplitude = wave;

