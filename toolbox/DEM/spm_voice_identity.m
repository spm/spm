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
% $Id: spm_voice_identity.m 7598 2019-05-25 13:09:47Z karl $


% Initial estimate of fundamental and formant frequencies
%--------------------------------------------------------------------------
global VOX
VOX.IT   = 1;
[Y,I,FS] = spm_voice_get_next(wfile);

if ~isfield(VOX,'F0')
    VOX.F0 = spm_voice_fundamental(Y,FS);
    VOX.F1 = 26 + VOX.F0/16;
end

%  estimate lexical and identity parameters under priors
%==========================================================================
SEG   = spm_voice_read(wfile,P,8);

%  accumulate evidence for reconsiders (assuming a Dirichlet distribution)
%--------------------------------------------------------------------------
L     = 0;
nw    = numel(SEG);
for w = 1:nw;
    L = L + SEG(w).L{3}/nw;
end

% posterior fundamental and formant  frequencies
%==========================================================================
ff0   = VOX.WHO(1).pE + VOX.R.F0;
ff1   = VOX.WHO(2).pE + VOX.R.F1;
F0    = exp(ff0'*L(:,1));
F1    = exp(ff1'*L(:,2));

% graphics if requested
%--------------------------------------------------------------------------
try VOX.formant;
    
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

