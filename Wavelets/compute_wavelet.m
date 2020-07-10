%
% computes the amplitud of the time-frequency decomposition using Morlet wavelets
%
% Peter Lakatos, modified by Julien Besle 02/12/2008

function [tf, frq, coi, W] = compute_wavelet(data,low_frq,high_frq,adrate, omega)

if ~ismember(nargin,[4 5]) 
    fprintf (1, 'Usage: [tf, frq, coi] = compute_wavelet(data, low_frq, high_frq, samp_rate,[omega]))\n');
    return;
end

% opengl 'neverselect';  % EDS2 3/27/2020 (Results in: "Warning: This option does not have any effect since it is no longer supported: neverselect ")
if ~exist('omega','var') || isempty(omega)      % the size of the wavelet, my correction only works with 6
    omega           = 6;                % a bigger value would give you more frequency
end                                     % accuracy but less time accuracy

num_tr = size(data, 2);   
WB = num_tr > 10;  % EDS2 3/27/2020 (waitbar not necessary for small sets)
if WB, waitb = waitbar(0); end  

if size(data,1) == 1 
    data = data';
end


for trial=num_tr:-1:1
	if WB  % EDS2 3/27/2020 (waitbar not necessary for small sets)
		a2=['trial ' num2str(num_tr-trial+1)];
		waitbar((num_tr-trial)/num_tr, waitb,a2);
	end
    [wave,period,scale,coi] = basewave4(data(:,trial),adrate,low_frq,high_frq,omega,0);            
    if nargout > 3, W(:, :, trial) = wave; end
    po=abs(wave);
    % the phase of the different frequency wavelets would be:
    % ph=angle(wave);
    if ~exist('tf', 'var')
        tf=zeros(size(po));
    end
    tf=tf+po;
end

tf=tf./trial;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HERE !!!!!
% 
% % Peter's correction (transformation gives a constant outcome when  varying the sample_rate (for a given number of samples), 
% tf=(tf/24)./ repmat(sqrt(period)',1,size(tf,2)) * sqrt(1000/adrate);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HERE !!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THERE !!!!
%
%my own correction so that 
%   1- the wavelet spectrum has the same scale as the the spectrum before multiplying by f
%   2- multiply by sqrt(f) so that the spectrum of a 1/f signal(pink noise) is flat
tf = tf*sqrt(5)*sqrt(size(tf,2))./repmat(sqrt(period)',1,size(tf,2));
%works fine for omega = 6 or 7
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THERE !!!!

frq=1./period;
coi = 1./coi;

if WB, close(waitb); end