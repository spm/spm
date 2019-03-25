function [I] = spm_voice_onsets(Y,FS)
% identifies intervals containing acoustic energy and post onset minima
% FORMAT [I] = spm_voice_onsets(Y,FS)
%
% Y    - timeseries
% FS   - sampling frequency
%
% i    - intervals (time bins) containing spectral energy
%
% This routine identifies epochs constaining spectral energy of the power
% envelope, defined as the law root mean square (RMS) power.the onset and
% offset of words is evaluated in terms of threshold crossings before and
% after the interquartile range (or the minimum if these do not exist).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onsets.m 7552 2019-03-25 10:46:03Z karl $

% find the interval that contains spectral energy
%==========================================================================

% identify threshold crossings in power
%--------------------------------------------------------------------------
n   = length(Y);
aY  = abs(Y);                                     % window absolute value
aY  = spm_conv(aY,FS/16);                         % smooth
aY  = aY - min(aY);                               % and normalise
aY  = aY/max(aY);
w   = ((1:n)'*2/n - 1).^2;
aY  = aY + w/8;

% find successive minima
%--------------------------------------------------------------------------
j   = find(diff(aY(1:end - 1)) < 0 & diff(aY(2:end)) > 0);
j   = j(j > FS/8);

% indices of interval containing spectral energy
%--------------------------------------------------------------------------
I     = {};
for i = 2:numel(j)
    k   = j(1):j(i);
    ni  = numel(k);
    if ni < 3*FS/4 && ni > FS/4;
        I{end + 1} = k;
    end
end


% graphics
%--------------------------------------------------------------------------
global VOX
if ~VOX.onsets; return, end

spm_figure('GetWin','onsets'); clf;

pst   = (1:n)/FS;
Ymax  = max(abs(Y));
subplot(2,1,1), plot(pst,Y), hold on
for i = 1:numel(j)
    plot([j(i),j(i)]/FS,[-1,1]*Ymax,'r')
end
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight, hold off
subplot(2,1,2), plot(pst,aY,'b'), hold on
for i = 1:numel(j)
    plot(pst(j(i)),aY(j(i)),'or')
end
title('Log energy','FontSize',16)
xlabel('peristimulus time (secs)'), spm_axis tight, hold off
drawnow

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










