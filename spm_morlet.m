function [C] = spm_morlet(s,k,wnum)
% Morlet wavelet transform (time-frequency analysis)
% FORMAT [C] = spm_morlet(s,k,wnum)
%
% s      - (t X n) time-series
% k      - Frequencies (cycles per time bin)
% wnum   - Wavelet number: default = 6
%
% C      - coefficents (complex)
%__________________________________________________________________________
%
% This routine returns a Morlet-like wavelet transform but uses a Hanning
% window, as opposed to a Gaussian window.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_morlet.m 6122 2014-07-25 13:48:47Z karl $


% setup and defaults
%--------------------------------------------------------------------------
if nargin < 3, wnum = 6;          end
if size(s,1) < size(s,2); s = s'; end

% window function (Hanning)
%--------------------------------------------------------------------------
[N,M] = size(s);
C     = zeros(N,length(k),M);

% s = randn(size(s))

% spectral density
%--------------------------------------------------------------------------
for i = 1:length(k)
    
    % Morlet convolution kernel
    %----------------------------------------------------------------------
    n     = 2*round(wnum/k(i)/2);
    h     = 1 - cos(2*pi*(1:n)/(n + 1));
    h     = h/sqrt(sum(h.^2));
    K     = h.*exp(-1j*(2*pi*k(i)*(0:(n - 1))));
    K     = spm_convmtx(K',N,'circular');

    % Transfer
    %----------------------------------------------------------------------
    C(:,i,:) = K*s;
    
end

return
