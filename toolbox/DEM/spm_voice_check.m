function [G,Y] = spm_voice_check(Y,FS,C)
% returns normalised spectral energy in acoustic range
% FORMAT [G,Y]  = spm_voice_check(Y,FS,C)
%
% Y    - timeseries
% FS   - sampling frequency
% C    - threshold [default: 1/16]
%
% Y    - high pass ( > 256 Hz) time series
% G    - root mean square energy of Y
%
% This routine identifies epochs constaining spectral energy in the
% acoustic range above of threshold for more than 200ms
% 
% see also spm_voice_filter.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_check.m 7574 2019-04-19 20:38:15Z karl $

% find the interval that contains spectral energy
%==========================================================================

% threshold - log ratio, relative to log(1) = 0;
%--------------------------------------------------------------------------
if nargin < 3, C = 1/16; end

% find periods of acoutic spectal energy > 200 ms
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















