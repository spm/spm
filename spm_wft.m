function [C] = spm_wft(s,k,n)
% Windowed fourier wavelet transform (time-frequency analysis)
% FORMAT [C] = spm_wft(s,k,n)
% s      - (t X n) time-series
% k      - Frequencies (cycles per window)
% n      - window length
% C      - (w X t X n) coefficients (complex)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2006-2022 Wellcome Centre for Human Neuroimaging


% window function (Hanning)
%--------------------------------------------------------------------------
[T,N] = size(s);
n     = round(n);
h     = 0.5*(1 - cos(2*pi*(1:n)/(n + 1)));
h     = h'/sum(h);
C     = zeros(length(k),T,N);


% spectral density
%--------------------------------------------------------------------------
for i = 1:length(k)
    W     = exp(-1j*(2*pi*k(i)*(0:(T - 1))/n))';
    for j = 1:N
        w        = conv(full(s(:,j)).*W,h);
        C(i,:,j) = w((1:T) + n/2);
    end
end
