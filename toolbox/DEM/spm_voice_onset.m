function [i] = spm_voice_onset(Y,FS,u,v)
% Identify intervals containing acoustic energy and post onset minima
% FORMAT [i] = spm_voice_onset(Y,FS,u,v)
%
% Y    - timeseries
% FS   - sampling frequency
% u,v  - thresholds for onset and offset [default: 1/16]
%
% i    - intervals (time bins) containing spectral energy
%
% This routine identifies epochs constaining spectral energy in the power
% envelope, defined as the root mean square (RMS) power. The onset and
% offset of words is evaluated in terms of the first and last threshold
% crossings.
%
% This routine is a simple version of spm_voice_onset and is retained for
% diagnostic purposes.
%
% see also: spm_voice_onsets.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onset.m 7750 2019-12-05 17:54:29Z spm $
% find the interval that contains spectral energy
%==========================================================================

% thresholds for onsets and offsets
%--------------------------------------------------------------------------
if nargin < 3, u = 1/16; v = 1/16; end

% identify threshold crossings in power
%--------------------------------------------------------------------------
n   = length(Y);                                  % length of time series
aY  = spm_hanning(n).*abs(Y);                     % window absolute value
aY  = spm_conv_full(aY,FS/16);                    % smooth
aY  = aY - min(aY);                               % and normalise
aY  = aY/max(aY);

i0  = find(aY > u,1,'first');                     % onset
iT  = find(aY > v,1,'last');                      % offsets

% indices of interval containing spectral energy
%--------------------------------------------------------------------------
i   = i0:iT;

% graphics
%--------------------------------------------------------------------------
global VOX
if ~VOX.onsets; return, end

spm_figure('GetWin','onsets'); clf;

pst   = (1:n)/FS;
Ymax  = max(abs(Y));
subplot(2,1,1), plot(pst,Y),      hold on
plot([i0,i0]/FS,[-1,1]*Ymax,'g'), hold on
plot([iT,iT]/FS,[-1,1]*Ymax,'r'), hold on
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight, hold off

subplot(2,1,2), plot(pst,aY,'b'),  hold on
plot(pst(iT),aY(iT),'or'),         hold on
plot(pst(i0),aY(i0),'og'),         hold on
plot(pst,u + spm_zeros(pst),':g'), hold on
plot(pst,v + spm_zeros(pst),':r'), hold off
title('Log energy','FontSize',16)
xlabel('peristimulus time (secs)'), spm_axis tight
drawnow

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)
