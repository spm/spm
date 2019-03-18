function [i,j] = spm_voice_onset(Y,FS,U)
% identifies intervals containing acoustic energy and post onset minima
% FORMAT [i,j] = spm_voice_onset(Y,FS,U)
%
% Y    - timeseries
% FS   - sampling frequency
% U    - threshold (log ratio) [default: 3]
%
% i    - intervals (time bins) containing spectral energy
% j    - indices of post onset minima
%
% This routine identifies epochs constaining spectral energy of the power
% envelope, defined as the law root mean square (RMS) power.the onset and
% offset of words is evaluated in terms of threshold crossings before and
% after the interquartile range (or the minimum if these do not exist).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onset.m 7546 2019-03-18 11:02:22Z karl $

% find the interval that contains spectral energy
%==========================================================================

% threshold - log ratio, relative to log(1) = 0;
%--------------------------------------------------------------------------
try
    u = 0 - U;
catch
    u = 0 - 4;
end

% fit fluctuations with a DCT and identify threshold crossings
%--------------------------------------------------------------------------
n     = length(Y);                           % length of time series
D     = spm_dctmtx(n,8);                     % discrete cosine transform 
fY    = D*(D'*log(abs(Y) + eps));            % fit log RMS
fY    = fY + hanning(n);                     % prior hanning 
fY    = fY - max(fY);                        % normalise to maximum
fY    = -6*fY/min(fY);                       % normalise range

% find interquartile range of spectral energy
%--------------------------------------------------------------------------
fY    = fY + hanning(n);                     % post hanning
sY    = cumsum(exp(fY));                     % cumulative RMS
w0    = find(sY > sY(n)/4,  1);              % first quartile
wT    = find(sY > sY(n)*3/4,1);              % third quartile

% identify threshold crossings
%--------------------------------------------------------------------------
i0    = find(fY(1:w0 - 1) < u,1,'last' );
iT    = find(fY(wT + 1:n) < u,1,'first');
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

% post onset minima if requested
%--------------------------------------------------------------------------
if nargout > 1
    j = find(diff(fY(1:end - 1)) < 0 & diff(fY(2:end)) > 0);
    j = j(j > (i0 + FS/4));
else
    j = [];
end

% graphics
%--------------------------------------------------------------------------
global VOX
if ~VOX.onsets; return, end

spm_figure('GetWin','onsets'); clf;

pst   = (1:n)/FS;
Ymax  = max(abs(Y));
subplot(2,1,1), plot(pst,Y),      hold on
plot([i0,i0]/FS,[-1,1]*Ymax,'r'), hold on
plot([iT,iT]/FS,[-1,1]*Ymax,'r'), hold on
plot([w0,w0]/FS,[-1,1]*Ymax,'g'), hold on
plot([wT,wT]/FS,[-1,1]*Ymax,'g'), hold off
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight

subplot(2,1,2), plot(pst,fY,'b'), hold on
plot(pst(j),fY(j),'or'),          hold on
plot(pst,u + spm_zeros(pst),':'), hold off
title('Log energy','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight
drawnow

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










