function [G] = spm_morlet_conv(G,k,wnum)
% temporal convolution of complex spectral responses with Morlet envelope
% FORMAT [G] = spm_morlet_conv(G,k,wnum)
%
% G      - (t x w x n x n) cross spectral density
% k      - Frequencies (cycles per time bin)
% wnum   - Wavelet number: default = 4
%
% G      - convolved cross spectral density
%__________________________________________________________________________
%
% This routine simply smooths a cross spectral response to emulate a 
% wavelet transform.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_morlet_conv.m 4768 2012-06-11 17:06:55Z karl $


% setup and defaults
%--------------------------------------------------------------------------
if nargin < 3, wnum = 4; end

% convolution
%--------------------------------------------------------------------------
for w = 1:size(G,2)
    n     = 2*round(wnum/k(w)/2);
    h     = 1 - cos(2*pi*(1:n)/(n + 1));
    h     = h/sqrt(sum(abs(h).^2));
    h     = h.^2';
    for i = 1:size(G,3)
        for j = 1:size(G,4)
            g          = conv(G(:,w,i,j),h);
            G(:,w,i,j) = g((1:size(G,1)) + n/2);
        end
    end
end
