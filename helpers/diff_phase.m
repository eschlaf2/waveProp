function dphi = diff_phase(phi, varargin)
% Angular differences (usage matches diff)
% dphi = diff_phase(phi, n=1, dim=1);

dphi = angle(exp(1j*diff(phi, varargin{:})));

end

%% Old

% if nargin < 2, dim = 1; end
% dim = dim - 1;
% phi = shiftdim(phi, dim);
% 
% if isreal(phi), phi = exp(1j*phi); end
% A = phi(1:end-2, :);
% B = phi(2:end-1, :);
% C = phi(3:end, :);
% dphi = angle(A.*conj(B) + B.*conj(C));
% dphi = shiftdim(dphi, ndims(phi)-dim);
