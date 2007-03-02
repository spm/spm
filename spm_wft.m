function [C] = spm_wft(s,k,n)
% Windowed fourier wavelet transform (time-frequency analysis)
% FORMAT [C] = spm_wft(s,k,n);
% s      - 1-D time-series
% k      - Frequencies (cycles per window)
% n      - window length
% C      - coefficents (complex)
%__________________________________________________________________________


% window function (Hanning)
%--------------------------------------------------------------------------
N     = length(s);
h     = 0.5*(1 - cos(2*pi*(1:n)/(n + 1)));
h     = h/sum(h);
C     = zeros(length(k),N);

% pad time-series
%--------------------------------------------------------------------------
s     = s(:)';

% spectral density
%-----------------------------------------------------------
for i = 1:length(k)
	W      = exp(-sqrt(-1)*(2*pi*k(i)*[0:(N - 1)]/n));
	w      = conv(full(s).*W,h);
	C(i,:) = w([1:N] + n/2);
end
