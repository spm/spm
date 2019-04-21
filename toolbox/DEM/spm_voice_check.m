function [G,Y] = spm_voice_check(Y,FS,C)
% returns normalised spectral energy in acoustic range
% FORMAT [G,Y]  = spm_voice_check(Y,FS,C)
%
% Y    - timeseries
% FS   - sampling frequency
% C    - standard deviation of spectral smoothing [default: 1/16 seconds]
%
% Y    - high pass ( > 512 Hz) time series
% G    - spectral envelope
%
% This routine applies a high pass filter by subtracting a smoothed version
% of the timeseries (to suppress frequencies of lesson 512 Hz). The
% absolute value of the resulting timeseriesis then convolved with a
% Gaussian kernel, specified by C.this returns the spectral envelope in
% terms of the root mean square energy (normalised to a minimum of zero).
% 
% see also: spm_voice_filter.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_check.m 7575 2019-04-21 16:47:39Z karl $

% find the interval that contains spectral energy
%==========================================================================

% standard deviation of spectral smoothing [default: 1/16 seconds]
%--------------------------------------------------------------------------
if nargin < 3, C = 1/16; end

% high pass filter and evaluate spectral envelope
%--------------------------------------------------------------------------
Y  = Y - spm_conv(Y,FS/512);
G  = spm_conv(abs(Y),FS*C);
G  = G - min(G);

return

% graphics
%--------------------------------------------------------------------------
pst = (1:numel(Y))/FS;
subplot(2,1,1), plot(pst,G)
title('Log energy','FontSize',16)
xlabel('peristimulus time'), spm_axis tight
drawnow















