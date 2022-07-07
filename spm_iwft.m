function [s] = spm_iwft(C,k,n)
% Inverse windowed Fourier transform - continuous synthesis
% FORMAT [s] = spm_iwft(C,k,n);
% s      - 1-D time-series
% k      - Frequencies (cycles per window)
% n      - window length
% C      - coefficients (complex)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging


% window function (Hanning)
%--------------------------------------------------------------------------
N     = size(C,2);
s     = zeros(1,N);
C     = conj(C);

% spectral density
%--------------------------------------------------------------------------
for i = 1:length(k)
    W      = exp(-sqrt(-1)*(2*pi*k(i)*[0:(N - 1)]/n));
    w      = W.*C(i,:);
    s      = s + real(w);
end
