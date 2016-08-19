function [G] = spm_morlet_conv(G,w,dt,wnum)
% temporal convolution of complex spectral responses with Morlet envelope
% FORMAT [G] = spm_morlet_conv(G,w,dt,wnum)
%
% G      - (t x w x n x n) cross spectral density
% w      - Frequencies (Hz)
% dt     - sampling interval (sec)
% wnum   - Wavelet number: default = 2  s.d. = wnum/(2*pi*w)
%
% G      - convolved cross spectral density
%__________________________________________________________________________
%
% This routine simply smooths a cross spectral response to emulate a 
% wavelet transform.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_morlet_conv.m 6857 2016-08-19 15:17:06Z karl $


% setup and defaults
%--------------------------------------------------------------------------
if nargin < 4, wnum = 8; end
[nt,nw,ni,nj] = size(G);
iw            = 1:nw;

% get (non-stationary) convolution matrix for frequencies
%--------------------------------------------------------------------------
f     = (-nw:nw)';
H     = zeros(nw,nw);
for i = 1:nw
    s      = w(i)/wnum;
    h      = exp(-f.^2/(2*s^2));
    h      = h(iw + nw - i);
    H(:,i) = h/sum(h);
end

% convolution over frequencies
%--------------------------------------------------------------------------
for i = 1:ni
    for j = 1:nj
        G(:,:,i,j) = G(:,:,i,j)*H;
    end
end

% convolution over time
%--------------------------------------------------------------------------
for k = 1:nw
    s     = wnum/(2*pi*w(k));
    t     = -(s*4):dt:(s*4);
    h     = exp(-t.^2/(2*s^2));
    h     = spm_convmtx(h',nt,'square');
    h     = diag(1./sum(h,2))*h;
    for i = 1:ni
        for j = 1:nj
            G(:,k,i,j) = h*G(:,k,i,j);
        end
    end
end

