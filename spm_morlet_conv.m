function [G] = spm_morlet_conv(G,w,dt,wnum)
% temporal convolution of complex spectral responses with Morlet envelope
% FORMAT [G] = spm_morlet_conv(G,w,dt,wnum)
%
% G      - (t x w x n x n) cross spectral density
% w      - Frequencies (Hz)
% dt     - sampling interval (sec)
% wnum   - Wavelet number: default = 2  s.d. = wnum * 1/w
%
% G      - convolved cross spectral density
%__________________________________________________________________________
%
% This routine simply smooths a cross spectral response to emulate a 
% wavelet transform.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_morlet_conv.m 5975 2014-05-07 18:07:42Z karl $


% setup and defaults
%--------------------------------------------------------------------------
if nargin < 4, wnum = 2; end

% get (non-stationary) convolution matrix for frequencies
%--------------------------------------------------------------------------
N     = length(w);
t     = (-N:N)'/(2*w(end));
ind   = 1:N;
for i = 1:N
    s      = wnum/w(i);
    h      = abs(fftshift(fft(exp(-t.^2/(2*s^2)))));
    h      = h(ind + N - i);
    H(:,i) = h/sum(h);
end

% convolution over frequencies
%--------------------------------------------------------------------------
for i = 1:size(G,3)
    for j = 1:size(G,4)
        G(:,:,i,j) = G(:,:,i,j)*H;
    end
end

% convolution over frequencies
%--------------------------------------------------------------------------
ind    = 1:size(G,1);
for k = 1:N
    s     = wnum/w(k);
    n     = s*4;
    t     = -n:dt:n;
    h     = exp(-t.^2/(2*s^2));
    h     = h/sum(h);
    for i = 1:size(G,3)
        for j = 1:size(G,4)
            g          = conv(G(:,k,i,j),h);
            G(:,k,i,j) = g(ind + round(length(t)/2));
        end
    end
end