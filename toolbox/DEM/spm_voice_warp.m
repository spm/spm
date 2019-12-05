function [I] = spm_voice_warp(Y,N)
% Resample a vector to normalise the phase at a particular frequency
% FORMAT [I] = spm_voice_warp(Y,N)
%
% Y    - timeseries
% N    - number of cycles (i.e., scale of normalisation)
%
% I    - resampling indices
%
% This auxiliary routine returns the indices of a vector that realigns the phase,
% following a Hilbert transform at a frequency of N cycles per vector
% length; i.e., warps the vector to normalise the phase at a specified
% scalable frequency
% 
% This routine is not actually used but is retained for reference
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_warp.m 7750 2019-12-05 17:54:29Z spm $


% find Sigma points (maxima of Hilbert transform)
%==========================================================================

% Fourier transform and frequencies
%--------------------------------------------------------------------------
fY    = fft(Y(:));
n     = length(fY);
w     = (1:n) - 1;
f0    = N;
s     = f0/4;

% bandpass filtered and rescale the phases
%--------------------------------------------------------------------------
bY    = fY.*exp(-(w(:) - f0).^2/(2*s^2));        % filter
sY    = real(ifft(bY));                          % filtered vector
I     = phase(hilbert(sY))/(2*pi);               % cumulative phase
I     = I*n/N;                                   % converted into bins

% check for overflow
%--------------------------------------------------------------------------
I     = I - I(1) + 1;
I     = round(max(min(I,n),1));

return

% graphics
%--------------------------------------------------------------------------
subplot(2,1,1);
plot((1:n),Y,':',(1:n),Y(I),'-',(1:n),sY,':',(1:n),sY(I),'--')
title('Original and warped vector'),xlabel('bins')

subplot(2,1,2)
plot((1:n),I,(1:n),(1:n),':')
title('Indices'),xlabel('bins')
