function [Y,W] = spm_voice_iff(xY)
% Inverse decomposition at fundamental frequency
% FORMAT [Y,W] = spm_voice_iff(xY)
%
% xY    -  cell array of word structures
% xY.W  -  parameters - lexical
% xY.P  -  parameters - prosody
% xY.R  -  parameters - speaker
% 
% xY.P.amp - log amplitude
% xY.P.dur - log duration (sec)
% xY.P.lat - log latency  (sec)
% xY.P.tim - timbre     (a.u.)
% xY.P.inf - inflection (a.u.)

% xY.R.F0  - fundamental frequency (Hz)
% xY.R.F1  - format frequency (Hz
%
% Y     - reconstructed timeseries
% W     - formants (time-frequency representation): Q = U*xY.W*V'
%      
%
% This routine recomposes a timeseries from temporal basis sets at the
% fundamental frequency. In other words, it applies the reverse sequence
% of inverse transforms implemented by spm_voice_ff.m
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_iff.m 7750 2019-12-05 17:54:29Z spm $

% defaults
%--------------------------------------------------------------------------
global VOX
if VOX.mute && ~nargout && ~VOX.graphics
    return
end
try, FS = VOX.FS; catch, FS  = spm_voice_FS; end    % sampling frequency
try, RF = VOX.RF; catch, RF  = 0;            end    % random fluctuations
try, rf = VOX.rf; catch, rf  = 0;            end    % random fluctuations


% recompose and play
%--------------------------------------------------------------------------
if numel(xY) > 1
    for s = 1:numel(xY)
        spm_voice_iff(xY(s));
    end
    return
end

% reconstitute sample points
%--------------------------------------------------------------------------
M  = exp(xY.P.amp);                          % amplitude
L  = exp(xY.P.lat);                          % latency (sec)
T  = exp(xY.P.dur);                          % duration (seconds)
S  = exp(xY.P.tim);                          % log timbre
G  = xY.P.pch;                               % log pitch
P  = xY.P.inf;                               % inflection

F0 = exp(xY.R.F0);                           % fundamental frequency (Hz)
F1 = exp(xY.R.F1);                           % formant frequency (Hz)

% reconstitute intervals
%--------------------------------------------------------------------------
ni = fix(T*F0);                              % number of intervals
D  = spm_dctmtx(ni,numel(P));                % basis set for inflection
dI = D*P(:)*sqrt(ni)/F0;                     % systemic fluctuations
dI = dI + rf*randn(ni,1)/F0;                 % random   fluctuations
I  = fix([1; FS*cumsum(dI)]);                % cumulative intervals

% reconstitute format coefficients
%--------------------------------------------------------------------------
Ni = 256;                                    % number of formant bins
nj = round(FS/F1);                           % interval length
W  = spm_voice_Q(xY.W,G,Ni,ni);              % log formant frequencies
W  = W + randn(Ni,ni)*RF;                    % add random fluctuations
Q  = exp(S*W);                               % formants
Q  = bsxfun(@rdivide,Q,hamming(ni)');        % unwindow
Q  = bsxfun(@minus,Q,min(Q));                % supress white noise

% reconstitute timeseries
%--------------------------------------------------------------------------
jj = 0:(2*nj);
Nj = numel(jj);
D  = spm_dctmtx(Nj,Ni*4);
D  = D*kron(speye(Ni,Ni),[1 0 -1 0]');
Y  = zeros(I(end) + 2*nj,1);
for j = 1:ni
    ii    = I(j)  + jj;
    Y(ii) = Y(ii) + D*Q(:,j);
end

% scale amplitude
%--------------------------------------------------------------------------
Y  = 4*M*Y/sum(std(Q));

% add latency
%--------------------------------------------------------------------------
Y  = [zeros(fix(L*FS),1); Y];

% play timeseries if requested
%--------------------------------------------------------------------------
if ~ VOX.mute && ~ nargout
    sound(Y,FS);
end

% graphics  if requested
%--------------------------------------------------------------------------
if ~ isfield(VOX,'graphics'), return, end

if VOX.graphics
    
    % figure
    %----------------------------------------------------------------------
    spm_figure('GetWin','Voice (graphics)'); clf;
    
    % peristimulus time (seconds) and plot
    %----------------------------------------------------------------------
    pst = (1:numel(Y))/FS;
    subplot(2,2,1), plot(pst,Y,[L,L],[-M M],':')
    xlabel('time (sec)'), ylabel('amplitude')
    title('Timeseries','FontSize',16), axis square, spm_axis tight
    
    subplot(2,2,2), imagesc((1:ni)/F0,1000*[-nj,nj]/FS,D*Q)
    axis square, xlabel('time (seconds)'), ylabel('time (ms)')
    title('Transients','FontSize',16), set(gca,'YLim',[-8 8])
    
    subplot(4,2,5), imagesc(xY.W), axis square
    xlabel('coefficients'), ylabel('coefficients')
    title('Parameters','FontSize',16)

    subplot(4,2,6), imagesc((1:ni)/F0,(1:Ni)*F1,Q)
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title('Spectral (log) energy','FontSize',16), drawnow

    subplot(4,2,7), plot(1./dI), axis square, spm_axis tight
    xlabel('time (intervals)'), ylabel('fundamental frequency')
    title('Inflection','FontSize',16)
    
    subplot(4,2,8), imagesc((1:ni)/F0,(1:Ni)*F1,W)
    xlabel('time (seconds)'), ylabel('Formants (Hz)')
    title('Spectral decomposition','FontSize',16), drawnow
    
end


return

% graphics for illustrations
%--------------------------------------------------------------------------
subplot(2,1,1), hold off
for j = 1:4
    ii    = I(j)  + jj;
    plot(ii,D*Q(:,j),'Color',spm_softmax(randn(3,1))), hold on
end
xlabel('time (bins)'), ylabel('transients')
title('Fundamental and first formant intervals','FontSize',16), drawnow
spm_axis tight, box off
