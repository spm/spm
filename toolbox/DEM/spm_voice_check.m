function [G,Y] = spm_voice_check(Y,FS,C)
% Return normalised spectral energy in acoustic range
% FORMAT [G,Y] = spm_voice_check(Y,FS,C)
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
% Gaussian kernel, specified by C. This returns the spectral envelope in
% terms of the root mean square energy (normalised to a minimum of zero).
% 
% see also: spm_voice_filter.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_check.m 7750 2019-12-05 17:54:29Z spm $


% find the interval that contains spectral energy
%==========================================================================

% standard deviation of spectral smoothing [default: 1/16 seconds]
%--------------------------------------------------------------------------
if nargin < 3, C = 1/16; end

% high pass filter and evaluate spectral envelope
%--------------------------------------------------------------------------
i = find(Y,1) + fix(FS*C);                 % deal with zero padding
i = min(i,numel(Y));                       % check for length
Y = Y - spm_conv_full(Y,FS/512);           % high pass filter
G = spm_conv_full(abs(Y),FS*C);            % root mean square power
G = G - min(G(i:end));                     % remove baseline power

return

% graphics
%--------------------------------------------------------------------------
pst = (1:numel(Y))/FS;
subplot(2,1,1), plot(pst,G)
title('Log energy','FontSize',16)
xlabel('peristimulus time'), spm_axis tight
drawnow
