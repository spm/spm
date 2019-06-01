function [F0,F1] = spm_voice_identity(wfile,P)
% Evaluates the fundamental and formant frequencies of a speaker
% FORMAT [F0,F1] = spm_voice_identity(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - prior probability of first word
%
% F0     - fundamental frequency
% F1     - expected format frequency
%
% This routine estimates the fundamental and formant frequencies based upon
% a spoken word source.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_identity.m 7600 2019-06-01 09:30:30Z karl $


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
        dt = (IS - VOX.IT)/spm_voice_FS;
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
SEG   = spm_voice_read(read(wfile),P,size(P,2));

% accumulate evidence for reconsiders (assuming a Dirichlet distribution)
%--------------------------------------------------------------------------
L     = 0;
nw    = numel(SEG);
for w = 1:nw;
    L = L + SEG(w).L{3}/nw;
end

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

if ~isfield(VOX,'formant'), return, end

% graphics if requested
%--------------------------------------------------------------------------
if VOX.formant
    
    % fundamental frequency
    %----------------------------------------------------------------------
    spm_figure('GetWin','Formants');
    
    FF0 = exp(ff0);
    FF1 = exp(ff1);
    
    subplot(2,2,1)
    plot(FF0,L(:,1),'or',FF0,L(:,1),':r',[F0 F0],[0 1],'b')
    xlabel('frequency (Hz)'), ylabel('likelihood')
    title(sprintf('Fundamental frequency (F0 = %0.0f)',F0),'FontSize',16)
    axis square, spm_axis tight, box off
    
    % first formant
    %----------------------------------------------------------------------
    
    subplot(2,2,2)
    plot(FF1,L(:,2),'or',FF1,L(:,2),':r',[F1 F1],[0 1],'b')
    xlabel('frequency (Hz)'), ylabel('likelihood')
    title(sprintf('1st formant frequency (F1 = %0.0f)',F1),'FontSize',16)
    axis square, spm_axis tight, box off, drawnow
    
end

