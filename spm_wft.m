function [C] = spm_wft(s,k,n)
% Windowed fourier wavelet transform (time-frequency analysis)
% FORMAT [C] = spm_wft(s,k,n)
% s      - (t X n) time-series
% k      - Frequencies (cycles per window)
% n      - window length
% C      - coefficents (complex)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_wft.m 5219 2013-01-29 17:07:07Z spm $


% transpose so that time runs down columns
%--------------------------------------------------------------------------
if size(s,1) < size(s,2); s = s'; end

% window function (Hanning)
%--------------------------------------------------------------------------
[N,M] = size(s);
h     = 0.5*(1 - cos(2*pi*(1:n)/(n + 1)));
h     = h'/sum(h);
C     = zeros(length(k),N,M);

% spectral density
%--------------------------------------------------------------------------
for i = 1:length(k)
    W      = exp(-1j*(2*pi*k(i)*(0:(N - 1))/n))';
    for j = 1:M
        w        = conv(full(s(:,j)).*W,h);
        C(i,:,j) = w((1:N) + n/2);
    end
end
