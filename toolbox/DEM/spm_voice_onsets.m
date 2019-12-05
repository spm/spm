function [I] = spm_voice_onsets(Y,FS,C,U)
% Identify intervals containing acoustic energy and post onset minima
% FORMAT [I] = spm_voice_onsets(Y,FS,C,U)
%
% Y    - timeseries
% FS   - sampling frequency
% C    - Convolution kernel [Default: 1/16 sec]
% U    - crossing threshold [Default: 1/8  a.u]

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
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_onsets.m 7750 2019-12-05 17:54:29Z spm $


% find the interval that contains spectral energy
%==========================================================================
global VOX

if nargin < 2, FS = VOX.FS; end                  % onset  threshold
if nargin < 4, U  = 1/8;    end                  % height threshold
if nargin < 3, C  = 1/16;   end                  % smoothing

% identify threshold crossings in power
%--------------------------------------------------------------------------
[G,Y] = spm_voice_check(Y,FS,C);                 % smooth envelope
G     = G/mean(G);                               % normalise to mean
n     = length(G);                               % number of bins

% find zero crossings (of U), for isolated words
%--------------------------------------------------------------------------
j0    = find(G(1:end - 1) < U & G(2:end) > U);
jT    = find(G(1:end - 1) > U & G(2:end) < U);

% find minima
%--------------------------------------------------------------------------
j     = find(diff(G(1:end - 1)) < 0 & diff(G(2:end)) > 0);
if numel(j)
    j = j([FS; diff(j)] > FS/32);                % remove minima < 1/32
    j = j(G(j) < 2);                             % remove minima > 2*mean
end

% find onsets preceded by silence
%--------------------------------------------------------------------------
k0    = [];
for k = 1:numel(j0)
    i = (j0(k) - FS/8):(j0(k) - 1);
    i = fix(i(i < n & i > 1));
    if all(G(i) < U) && j0(k) < FS/2
        k0(end + 1) = j0(k);
    end
end

% find the last onset preceded by silence
%--------------------------------------------------------------------------
if isempty(k0), j0 = []; else, j0 = k0(end); end

% find offsets followed by silence
%--------------------------------------------------------------------------
kT    = [];
for k = 1:numel(jT)
    i = (jT(k) + 1):(jT(k) + FS/8);
    i = fix(i(i < n & i > 1));
    if all(G(i) < U) && jT(k) > FS/2
        kT(end + 1) = jT(k);
    end
end

% find the first offset followed by silence
%--------------------------------------------------------------------------
if isempty(kT), jT = []; else, jT = kT(1); end


% use the first minima in the absence of zero crossings
%--------------------------------------------------------------------------
if isempty(j0)
    try
        j0 = j(1);
    catch
        j0 = 1;
    end
end

% use the last minima in the absence of zero crossings
%--------------------------------------------------------------------------
if isempty(jT)
    try
        jT = j(end);
    catch
        jT = n;
    end
end

% use the last sample if offset precedes onset
%--------------------------------------------------------------------------
if (jT - j0) < 1, jT = n; end

% add internal minima
%--------------------------------------------------------------------------   
i   = j(j < jT(end) & j > (j0 + FS/16));
jT  = sort(unique([jT; i]));


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
        if ni > FS/16 && ni < FS
            I{end + 1} = k;
        end
    end
end

% sort lengths (longest last), with 3 minima or less
%--------------------------------------------------------------------------
for i = 1:numel(I)
    ni(i) = numel(I{i});
end
[d,j] = sort(ni,'ascend');
I     = I(j(1:min(3,end)));

% graphics(if requested)
%==========================================================================
if ~VOX.onsets
    return
else
    spm_figure('GetWin','onsets'); clf;
end

% timeseries
%--------------------------------------------------------------------------
pst   = (1:n)/FS;
subplot(2,1,1)
plot(pst,Y/max(Y), 'b'),     hold on
plot(pst,G/max(G),':b'),     hold on
plot([1 1]/2,[-1 1],'b'),    hold on
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
for i = 1:numel(j0), plot(pst(j0(i)),G(j0(i)),'og'), end
for i = 1:numel(jT), plot(pst(jT(i)),G(jT(i)),'or'), end
title('Spectral envelope','FontSize',16)
xlabel('peristimulus time (secs)'), spm_axis tight, hold off
drawnow, pause(1/4)

% uncomment to play interval
%--------------------------------------------------------------------------
% sound(Y(i),FS)
