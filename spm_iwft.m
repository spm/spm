function [s] = spm_iwft(C,k,n)
% Inverse windowed Fourier transform - continuous synthesis
% FORMAT [s] = spm_iwft(C,k,n);
% s      - 1-D time-series
% k      - Frequencies (cycles per window)
% n      - window length
% C      - coefficents (complex)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_iwft.m 1131 2008-02-06 11:17:09Z spm $


% window function (Hanning)
%--------------------------------------------------------------------------
N     = size(C,2);
s     = zeros(1,N);
C     = conj(C);

% spectral density
%-----------------------------------------------------------
for i = 1:length(k)
    W      = exp(-sqrt(-1)*(2*pi*k(i)*[0:(N - 1)]/n));
    w      = W.*C(i,:);
    s      = s + real(w);
end

