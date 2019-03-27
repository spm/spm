function [I] = spm_voice_onsets(Y,FS,U,u)
% identifies intervals containing acoustic energy and post onset minima
% FORMAT [I] = spm_voice_onsets(Y,FS,U,u)
%
% Y    - timeseries
% FS   - sampling frequency
% U    - onset threshold [Defult: 1/16]
% u    - crossing threshold [Defult: 1/16]
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
% $Id: spm_voice_onsets.m 7557 2019-03-27 17:11:16Z karl $

% find the interval that contains spectral energy
%==========================================================================
if nargin < 3, U = 1/16; end
if nargin < 4, u = 1/16; end

% identify threshold crossings in power
%--------------------------------------------------------------------------
aY  = spm_conv(abs(Y),FS/16);                         % smooth
n   = length(aY);
aY  = aY - min(aY);
aY  = aY/max(aY);
y0  = aY;
yT  = aY;

% find successive minima
%--------------------------------------------------------------------------
j0    = find(y0(1:end - 1) < u & y0(2:end) > u);
jT    = find(yT(1:end - 1) > u & yT(2:end) < u);

% boundary conditions
%--------------------------------------------------------------------------
if numel(j0) > 1
    j0 = j0(j0 > FS*U & j0 < FS/2);
end
if numel(j0) > 1
    j0 = j0(1);
end
if numel(jT) > 1
    jT = jT(jT > FS/2);
end
if numel(jT) < 1
    jT = n;
end

% indices of interval containing spectral energy
%--------------------------------------------------------------------------
I     = {};
for i = 1:numel(j0)
    for j = 1:numel(jT)
        k   = j0(i):jT(j);
        ni  = numel(k);
        if ni < 3*FS/4 && ni > FS/8;
            I{end + 1} = k;
        end
    end
end

% sort lengths
%--------------------------------------------------------------------------
for i = 1:numel(I)
    ni(i) = numel(I{i});
end
[i,j] = sort(ni,'descend');
I     = I(j);

% graphics
%--------------------------------------------------------------------------
global VOX
if ~VOX.onsets; return, end

spm_figure('GetWin','onsets'); clf;

pst   = (1:n)/FS;
Y     = Y/max(Y);
subplot(2,1,1)
plot(pst,Y,'b'),             hold on
plot(pst,aY/max(aY),':b'),   hold on
plot(pst,Y,'b'),             hold on
plot(pst,aY,':b'),           hold on
plot(pst,0*pst + u,'-.'),    hold on
plot([1 1]/2,[-1 1],'b'),    hold on
plot([1 1]*U,[-1 1],'--g'),  hold on
for i = 1:numel(I)
    x = [I{i}(1),I{i}(end),I{i}(end),I{i}(1)]/FS;
    y = [-1,-1,1,1];
    c = spm_softmax(rand(3,1))';
    h = fill(x,y,c);
    set(h,'Facealpha',1/8,'EdgeAlpha',1/8);
end
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight, hold off

subplot(2,1,2)
plot(pst,aY,'b'),            hold on
plot(pst,0*pst + u,'-.'),    hold on
plot([1 1]/2,[0 1],'b'),     hold on
plot([1 1]*U,[0 1],'--g'),   hold on
for i = 1:numel(j0), plot(pst(j0(i)),y0(j0(i)),'og'), end
for i = 1:numel(jT), plot(pst(jT(i)),yT(jT(i)),'or'), end
title('Log energy','FontSize',16)
xlabel('peristimulus time (secs)'), spm_axis tight, hold off
drawnow, pause(1/4)

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










