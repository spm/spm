function [G] = spm_voice_check(Y,FS,C)
% returns normalised spectral energy in acoustic range
% FORMAT [G] = spm_voice_check(Y,FS,C)
%
% Y    - timeseries
% FS   - sampling frequency
% C    - threshold [default: 1/16]
%
% i    - intervals (time bins) containing spectral energy
%
% This routine identifies epochs constaining spectral energy in the
% acoustic range above of threshold for more than 200ms
% 
% see also spm_voice_filter.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_check.m 7566 2019-04-03 12:15:50Z karl $

% find the interval that contains spectral energy
%==========================================================================

% threshold - log ratio, relative to log(1) = 0;
%--------------------------------------------------------------------------
if nargin < 3, C = 1/16; end

% find periods of acoutic spectal energy > 200 ms
%--------------------------------------------------------------------------
Y  = Y - spm_conv(Y,FS/256);
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















