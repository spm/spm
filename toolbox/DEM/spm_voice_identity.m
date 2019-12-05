function [F0,F1] = spm_voice_identity(wfile,P)
% Evaluate the fundamental and formant frequencies of a speaker
% FORMAT [F0,F1] = spm_voice_identity(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - prior probability of first word
%
% F0     - fundamental frequency
% F1     - expected format frequency
%
% NB: automatically updates VOX.F0 and VOX.F1 when called
%
% This routine estimates the fundamental and formant frequencies based upon
% a spoken word source. This routine is used in conjunction with
% spm_voice_fundamental to provide a more refined estimate of fundamental
% and first formant frequencies based upon speech with known lexical
% content.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_identity.m 7750 2019-12-05 17:54:29Z spm $


% global VOX
%--------------------------------------------------------------------------
global VOX

% get source (recorder) and FS
%--------------------------------------------------------------------------
[FS,read] = spm_voice_FS(wfile);

% ensure 2 second of data has been accumulated
%--------------------------------------------------------------------------
if isa(wfile,'audiorecorder')
    IS = get(wfile,'TotalSamples');
    if ~IS
        stop(VOX.audio);
        record(VOX.audio,8);
        pause(4);
    else
        dt = (IS - VOX.IT)/FS;
        pause(4 - dt);
    end
end

% Initial estimate of fundamental and formant frequencies
%--------------------------------------------------------------------------
if ~isfield(VOX,'F0')
    VOX.F0 = spm_voice_fundamental(read(wfile),FS);
    VOX.F1 = 26 + VOX.F0/16;
end

% estimate lexical and identity parameters under priors
%==========================================================================
SEG   = spm_voice_read(read(wfile),P);

% accumulate evidence for reconsiders (assuming a Dirichlet distribution)
%--------------------------------------------------------------------------
L     = 0;
nw    = numel(SEG);
for w = 1:nw
    L = L + SEG(w).L{3};
end
L     = spm_softmax(L/8);

% posterior fundamental frequencies
%==========================================================================
ff0   = VOX.WHO(1).pE + VOX.R.F0;
ff1   = VOX.WHO(2).pE + VOX.R.F1;
F0    = exp(ff0'*L(:,1));

% with epirical prior for the first formant
%--------------------------------------------------------------------------
Ef    = log(26 + F0/16);                      % prior expectation
Pf    = 256;                                  % prior precision
PL    = -(ff1 - Ef).^2*Pf;                    % prior potential
PL    = spm_softmax(PL + log(L(:,2)));        % posterior probability
F1    = exp(ff1'*PL);                         % posterior expectation

% update VOX structure
%--------------------------------------------------------------------------
VOX.F0 = F0;
VOX.F1 = F1;

if ~isfield(VOX,'formant'), return, end

% graphics if requested
%--------------------------------------------------------------------------
if VOX.formant
    
    % fundamental frequency
    %----------------------------------------------------------------------
    spm_figure('GetWin','Frequencies');
    
    FF0 = exp(ff0);
    FF1 = exp(ff1);
    
    subplot(3,2,5)
    plot(FF0,L(:,1),'og',FF0,L(:,1),':r',[F0 F0],[0 1],'b')
    xlabel('frequency (Hz)'), ylabel('likelihood')
    title(sprintf('Fundamental (F0 = %0.0f Hz)',F0),'FontSize',16)
    axis square, spm_axis tight, box off
    
    % first formant
    %----------------------------------------------------------------------
    subplot(3,2,6)
    plot(FF1,L(:,2),'og',FF1,L(:,2),':r',[F1 F1],[0 1],'b')
    xlabel('frequency (Hz)'), ylabel('likelihood')
    title(sprintf('1st formant  (F1 = %0.0f (%0.0f) Hz)',F1,F1*32),'FontSize',16)
    axis square, spm_axis tight, box off, drawnow
    
    % supplement spectral plots
    %----------------------------------------------------------------------
    subplot(3,1,1), hold on
    a  = get(gca,'YLim'); a = a(2);
    plot([F0 F0],[0 a],'g')
    for j = 1:4
        plot([F0 F0]*j,[0 a],'-.g')
    end, hold off
    
    subplot(3,1,2), hold on
    for j = 1:8
        plot(16*[F1 F1]*j,[0 a],'g')
        for k = 1:8
            plot(16*[F1 F1]*(j - 1 + (k - 1)/8),[0 a],':g')
        end
    end, hold off
    
end
