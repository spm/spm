function [C] = spm_mwt(s,k,w)
% Morlet wavelet transform (time-frequency analysis)
% FORMAT [C] = spm_mwt(s,k,w)
% s      - (T X N) time-series
% k      - Frequencies (cycles per window)
% w      - window length
% C      - (w X T X N) coefficients (complex)
%
% see also: spm_iwft
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2006-2022 Wellcome Centre for Human Neuroimaging


% window function (Hanning)
%--------------------------------------------------------------------------
[T,N] = size(s);


% spectral density
%--------------------------------------------------------------------------
C     = zeros(length(k),T,N);
for i = 1:length(k)
    n     = w*5/k(i);
    n     = round(min(n,w*32));
    h     = 0.5*(1 - cos(2*pi*(1:n)/(n + 1)));
    h     = h'/sum(h + eps);
    W     = exp(-1j*(2*pi*k(i)*(0:(T - 1))/n))';
    for j = 1:N
        C(i,:,j) = conv(full(s(:,j)).*W,h,'same');
    end
end

return
