function [i] = spm_voice_onset(Y,FS)
% decomposition at fundamental frequency
% FORMAT [i] = spm_voice_onset(Y,FS)
%
% Y    - timeseries
% FS   - sampling frequency
% F0   - fundamental frequency (glottal pulse rate)
%
% I    - intervals (time bins): mean(I) = DI = FS/F0
% DJ   - FS/F1: F1   - formant frequency
%
% xY.Y  - timeseries
% xY.y  - unwrapped timeseries
% xY.FS - sampling frequency (hertz)
% xY.ni - fundamental frequency(hertz)
% xY.ti - duration (seconds)
% xY.ci - coefficients
%
% This routine decomposes a timeseries into a temporal basis set at the
% fundamental frequency
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onset.m 7528 2019-02-06 19:19:49Z karl $

% find the interval that contains spectral energy
%==========================================================================

% find interquartile range of spectral energy
%--------------------------------------------------------------------------
n     = length(Y);                           % Length of time series
aY    = abs(Y);                              % root power
sY    = cumsum(hanning(n).*sqrt(aY));        % cumulative (root) power
w0    = find(sY > sY(n)/4,  1);              % first quartile
wT    = find(sY > sY(n)*3/4,1);              % third quartile
aY    = spm_conv(aY,FS/16);                  % smooth and normalise
aY    = aY - min(aY);
Ymax  = max(aY);
aY    = aY/Ymax;

% fit fluctuations with a DCT and identify threshold crossings
%--------------------------------------------------------------------------
D     = spm_dctmtx(n,8);                     % discrete cosine transform 
lY    = log(aY + exp(-4));                   % log transformed
fY    = D*(D'*lY);                           % fit
u     = -3;                                  % threshold (nats)
i0    = find((fY(1:w0 - 1) < u) & (fY(2:w0)     > u),1,'last' );
iT    = find((fY(wT:n - 1) > u) & (fY(wT + 1:n) < u),1,'first');
iT    = iT + wT;

% use the minima before and after interquartile range, if necessary
%--------------------------------------------------------------------------
if isempty(i0)
    [d,i0] = min(fY(1:w0));
end
if isempty(iT)
    [d,iT] = min(fY(wT:n));
    iT     = iT + wT - 1;
end

% indices of interval containing spectral energy
%--------------------------------------------------------------------------
i     = i0:iT;

global voice_options
if ~voice_options.onsets;
    return
end

% graphics
%--------------------------------------------------------------------------
pst   = (1:n)/FS;
Ymax  = max(abs(Y));
subplot(2,1,1), plot(pst,Y),      hold on
plot([i0,i0]/FS,[-1,1]*Ymax,'r'), hold on
plot([iT,iT]/FS,[-1,1]*Ymax,'r'), hold on
plot([w0,w0]/FS,[-1,1]*Ymax,'g'), hold on
plot([wT,wT]/FS,[-1,1]*Ymax,'g'), hold off
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight

subplot(2,1,2), plot(pst,lY,'c',pst,fY,'b'), hold on
plot(pst,u + spm_zeros(pst),':'), hold off
title('Log energy','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










