function [i] = spm_voice_check(Y,FS,U)
% decomposition at fundamental frequency
% FORMAT [i] = spm_voice_onset(Y,FS,U)
%
% Y    - timeseries
% FS   - sampling frequency
% U    - threshold [default: 3]
%
% i    - intervals (time bins) containing spectral energy
%
% This routine identifies epochs constaining spectral energy in the
% acoustic range above of threshold for more than 200ms
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_check.m 7535 2019-03-03 20:27:09Z karl $

% find the interval that contains spectral energy
%==========================================================================

% threshold - log ratio, relative to log(1) = 0;
%--------------------------------------------------------------------------
if nargin < 3, U = 3; end

% find periods of acoutic spectal energy > 200 ms
%--------------------------------------------------------------------------
G  = spm_voice_filter(Y,FS);
G  = log(G);
G  = G - min(G);
if sum(G > U) < FS/5
    i  = [];
    return
else
    i = spm_voice_onset(Y,FS);
end

return

% graphics
%--------------------------------------------------------------------------
subplot(2,1,1), plot(G),    hold on
plot(U + spm_zeros(G),':'), hold off
title('Log energy','FontSize',16)
xlabel('peristimulus time'), spm_axis tight
drawnow










