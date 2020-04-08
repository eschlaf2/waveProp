% BASEWAVE2 convolves a Morlet-based wavelet with the data to compute the
% time-frequency representation of any signal.  The function returns the
% complex number representation of these signals, and it can be utilized as
% follows:
%
%   [WAVE,PERIOD,SCALE,COI] = BASEWAVE2(DATA,Fs,LOWF,HIF,WAITH), 
%
% where the variables are:
%
%   WAVE - time-frequency representation of the signal (power and angle can
%          be calculated from these complex numbers);
%   PERIOD - inverse of the frequency scale;
%   SCALE - ?
%   COI - cone of influence indicates where the wavelet analysis is skewed
%         because of edge artifacts;
%   DATA - the signal in the time domain;
%   adrate - sampling frequency;
%   Laxis - lower frequency of the range for which the transform is to be
%          done;
%   Naxis - higher frequency of the range for which the transform is to be
%         done;
%   WAITH - a handle to the qiqi waitbar.
%
% See also BASEWAVE, ABRA, WAVELET.

% Original code written by Torrence
% Modified by Peter Lakatos many times
% Modified by Ankoor Shah to return the complex representation


function [wave,period,scale,coi]=basewave4(Y,adrate,Laxis,Naxis,k0,waitc);

dt      = 1/adrate;		   
dj      = 0.08;
s0      = 1/(Naxis+(0.1*Naxis));    % the smallest resolvable scale
pad     = 0;                        % zero padding (1-on, 0-off)
% k0      = 6;                        % the mother wavelet parameter (wavenumber), default is 6.
frq_lo  = Laxis;
frq_hi  = Naxis;
n = length(Y);
J1 = (log2(n*dt/s0))/dj;                  % [Eqn(10)] (J1 determines the largest scale)

%....waitbar
if waitc==0;
else waitbar(0, waitc);
end

%....az adatsor atalakitasa
n1 = length(Y); 
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];                               % k n/2 elembol all
k = k*((2*pi)/(n*dt));                  
k = [0, k, -k(fix((n-1)/2):-1:1)];              % k ismet n elembol all - angular frequency


%....compute FFT of the (padded) time series (DFT)
f = fft(x);                                     % [Eqn(3)]

%....construct SCALE array

scale = s0*2.^((0:J1)*dj);                      % [Eqn(9)]  (choice of scales to use in wavelet transform ([Eqn(9)]))
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));  % Scale-->Fourier [Sec.3h]
coi = fourier_factor/sqrt(2);                   % Cone-of-influence [Sec.3g]
period = fourier_factor*scale;
xxx=min(find(1./period<frq_lo));
period=fliplr(period(1:xxx));
scale=fliplr(scale(1:xxx));
%....find freqency-equivalent scales
% aa=1;
% for a=frq_lo:0.2:frq_hi
%     ii(aa)=max(find(period<(1/a)));
%     scale1(aa)=scale(ii(aa));
%     aa=aa+1;
% end
% 
% fscale=aa-1;

scale1=scale;
fscale=size(period,2);

%....construct empty WAVE array
wave = zeros(fscale,n);     % define the wavelet array
wave = wave + i*wave;       % make it complex

%....loop through all scales and compute transform
for a1 = 1:fscale
            if waitc==0;
            else waitbar(a1/fscale, waitc);
            end
    expnt = -(scale1(a1).*k - k0).^2/2.*(k > 0.);
    norm = sqrt(scale1(a1)*k(2))*(pi^(-0.25))*sqrt(n);      % total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt);
    daughter = daughter.*(k > 0.);                          % Heaviside step function
	wave(a1,:) = ifft(f.*daughter);                         % wavelet transform[Eqn(4)]
    clear expnt,daughter;
end

period = fourier_factor*scale1;                                 % az �jfajta, fourier frekvenci�knak megfelelo period
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];   % COI [Sec.3g]

% wave = wave(Laxis:fscale,1:n1);                                  % get rid of padding before returning
% period = period(Laxis:fscale);

%powerx = (abs(wave)).^2;                                        % compute wavelet power spectrum
%powerx = (powerx*1000)/adrate;                                  % ez egy korrekci�, ami a torrence-ben nem volt benne

return