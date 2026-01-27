function [s] = spm_imwt(C,k,w)
% Inverse windowed Fourier transform - continuous synthesis
% FORMAT [s] = spm_imwt(C,k,w);
% s      - 1-D time-series
% k      - Frequencies (cycles per window)
% w      - window length
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
    n = w*5/k(i);
    n = round(min(n,w*32));
    W = exp(-sqrt(-1)*(2*pi*k(i)*(0:(N - 1))/n));
    s = s + real(W.*C(i,:));
end
