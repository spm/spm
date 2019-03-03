function [i,j] = spm_voice_onset(Y,FS,U)
% decomposition at fundamental frequency
% FORMAT [i,j] = spm_voice_onset(Y,FS,U)
%
% Y    - timeseries
% FS   - sampling frequency
% U    - threshold (log ratio) [default: 3]
%
% i    - intervals (time bins) containing spectral energy
% j    - indices of post onset minima
%
% This routine identifies epochs constaining spectral energy
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onset.m 7535 2019-03-03 20:27:09Z karl $

% find the interval that contains spectral energy
%==========================================================================

% threshold - log ratio, relative to log(1) = 0;
%--------------------------------------------------------------------------
try
    u = 0 - U;
catch
    u = 0 - 3;
end
    
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
i     = i0:iT;                                       % interval

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
global voice_options
if ~voice_options.onsets; return, end

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
plot(pst(j),fY(j),'or'),          hold on
plot(pst,u + spm_zeros(pst),':'), hold off
title('Log energy','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight
drawnow

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










