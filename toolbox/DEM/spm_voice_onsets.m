function [I] = spm_voice_onsets(Y,FS,T,U,C)
% identifies intervals containing acoustic energy and post onset minima
% FORMAT [I] = spm_voice_onsets(Y,FS,T,U,C)
%
% Y    - timeseries
% FS   - sampling frequency
% T    - onset    threshold [Default: 1/16 sec]
% U    - crossing threshold [Default: 1/16 a.u]
% C    - Convolution kernel [Default: 1/16 sec]
%
% I{i} - cell array of intervals (time bins) containing spectral energy
%
% This routine identifies epochs constaining spectral energy of the power
% envelope, defined as the root mean square (RMS) power. The onset and
% offset of words is evaluated in terms of threshold crossings before and
% after the midpoint of a one second epoch. These are supplemented with
% internal minima (after the spectral peak).
%
% see also: spm_voice_onset.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onsets.m 7583 2019-05-02 12:10:08Z karl $

% find the interval that contains spectral energy
%==========================================================================
global VOX

if nargin < 2, FS = VOX.FS; end                       % onset  threshold
if nargin < 3, T  = 1/8;    end                       % onset  threshold
if nargin < 4, U  = 1/1024; end                       % height threshold
if nargin < 5, C  = 1/16;   end                       % smoothing

% identify threshold crossings in power
%--------------------------------------------------------------------------
[G,Y] = spm_voice_check(Y,FS,C);                      % smooth RMS

% find zero crossings (of U)
%--------------------------------------------------------------------------
j0    = find(G(1:end - 1) < U & G(2:end) > U,1);
jT    = find(G(1:end - 1) > U & G(2:end) < U);

% boundary conditions
%--------------------------------------------------------------------------
j0    = j0(j0 > FS*T & j0 < FS/2);
jT    = jT(jT > FS/2);

if numel(j0) < 1,  j0 = fix(FS*T); end
if numel(jT) < 1,  jT = length(G); end


% and supplement offsets with (internal) minima
%--------------------------------------------------------------------------
i     = find(diff(G(1:end - 1)) < 0 & diff(G(2:end)) > 0);
i     = i(G(i) < max(G)/2 & i < jT(end) & i > FS/2);
jT    = sort(unique([jT; i]));

% check for silence after offset (and remove internal minima)
%--------------------------------------------------------------------------


% indices of interval containing spectral energy
%--------------------------------------------------------------------------
I     = {};
for i = 1:numel(j0)
    for j = 1:numel(jT)
        
        % k-th interval
        %------------------------------------------------------------------
        k  = j0(i):jT(j);
        ni = numel(k);
        
        % retain intervals of plausible length
        %------------------------------------------------------------------
        if ni > FS/16
            I{end + 1} = k;
        end
    end
end

% sort lengths (longest first)
%--------------------------------------------------------------------------
for i = 1:numel(I)
    ni(i) = numel(I{i});
end
[d,j] = sort(ni,'descend');
I     = I(j);

% graphics(if requested)
%==========================================================================
if ~VOX.onsets
    return
else
    spm_figure('GetWin','onsets'); clf;
end

% timeseries
%--------------------------------------------------------------------------
pst   = (1:length(G))/FS;
subplot(2,1,1)
plot(pst,Y/max(Y), 'b'),     hold on
plot(pst,G/max(G),':b'),     hold on
plot([1 1]/2,[-1 1],'b'),    hold on
plot([1 1]*T,[-1 1],'--g'),  hold on
for i = 1:numel(I)
    x = [I{i}(1),I{i}(end),I{i}(end),I{i}(1)]/FS;
    y = [-1,-1,1,1];
    c = spm_softmax(rand(3,1))';
    h = fill(x,y,c);
    set(h,'Facealpha',1/8,'EdgeAlpha',1/8);
end
title('Onsets and offsets','FontSize',16)
xlabel('peristimulus time (seconds)'), spm_axis tight, hold off

% envelope and threshold crossings
%--------------------------------------------------------------------------
subplot(2,1,2)
plot(pst,G,'r'), hold on
plot(pst,0*pst + U,'-.'),         hold on
plot([1 1]/2,[0 max(G)],'b'),     hold on
plot([1 1]*T,[0 max(G)],'--g'),   hold on
for i = 1:numel(j0), plot(pst(j0(i)),G(j0(i)),'og'), end
for i = 1:numel(jT), plot(pst(jT(i)),G(jT(i)),'or'), end
title('Spectral envelope','FontSize',16)
xlabel('peristimulus time (secs)'), spm_axis tight, hold off
drawnow, pause(1/4)

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)










