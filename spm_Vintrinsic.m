function [V] = spm_Vintrinsic(n,RT,form,Vpar)
% Computes an assumed intrinsic autocorelation function
% FORMAT [V] = spm_Vintrinsic(n,form,Vparam)
%
% n      - sessions number of scans
% RT	 - Repetition time in seconds
% form	 - '1/f' 
% Vpar	 - Parameters pertaining to the form
%
% V	 - (intrinsic) autoccorrelation matrix
%__________________________________________________________________________
% %W% Jean-Baptiste Poline, Karl Friston %E%

%--------------------------------------------------------------------------
if nargin < 4 error('not enough arguments in spm_Vintrinsic'); end;


% frequencies
%--------------------------------------------------------------------------
if rem(n,2) == 0
	dF  = 1/(2*n*RT);
else 
	dF  = 1/(2*(n - 1)*RT);
end
F    = [1:n]*dF;


% Spectral density assumed
%--------------------------------------------------------------------------
switch form

	case('1/f')
	%------------------------------------------------------------------
	G   = (Vpar(1)./F + 1).^2;

	otherwise, error('unrecognised type')
end

% autocorrelation functional and matrix
%--------------------------------------------------------------------------
R    = fftshift(real(ifft(G([1:n n:-1:2]))));
R    = R/max(R);
V    = toeplitz(R([1:n] + n - 1));

