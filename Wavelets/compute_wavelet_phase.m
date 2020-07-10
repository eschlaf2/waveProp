%
% computes the phase of the time-frequency decomposition using Morlet wavelets
%
% Peter Lakatos, modified by Julien Besle 02/12/2008

function [tf, frq, coi] = compute_wavelet_phase(data,low_frq,high_frq,adrate)

if nargin ~= 4 | nargout~=3 
    fprintf (1, 'Usage: [tf, frq, coi] = compute_wavelet_phase(data, low_frq, high_frq, samp_rate))\n');
    return;
end

opengl 'neverselect';
omega           = 6;        % the size of the wavelet, my correction only works with 6
                            % a bigger value would give you more frequency
                            % accuracy but less time accuracy

waitb = waitbar(0);

if size(data,1) == 1 
    data = data';
end

for trial=1:size(data,2)
    a2=['trial ' num2str(trial)];
    waitbar(trial/size(data,2), waitb,a2);
    [wave,period,scale,coi] = basewave4(data(:,trial),adrate,low_frq,high_frq,omega,0);            
    
    ph=angle(wave);
    
    if trial==1
        tf=zeros(size(ph));
    end
    tf=tf+ph;
end

tf=tf./trial;
frq=1./period;
coi = 1./coi;
close(waitb)