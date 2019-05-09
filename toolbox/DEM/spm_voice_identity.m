function [L,F0,F1] = spm_voice_identity(wfile,P)
% Evaluates the fundamental and formant frequencies of a word
% FORMAT [L,F0,F1] = spm_voice_identity(wfile,P)
%
% wfile  - .wav file, audiorecorder object or (double) time series
% P      - prior probability of first word
%
% L      - log-likelihood over F1 (identity attribute)
% F0     - fundamental frequency
% F1     - expected format frequency
%
% requires the following in the global variable VOX:
%
% LEX    - lexical structure array
% PRO    - prodidy structure array
% WHO    - speaker structure array
% FS     - sampling    frequency (Hz)
% F0     - fundamental frequency (Hz)
% IT     - index or pointer to offset of last word (i.e., CurrentSample)
%
% This routine estimates the fundamental and formant frequencies based upon
% a spoken word source.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_identity.m 7589 2019-05-09 12:57:23Z karl $


%% log prior over lexical content
%==========================================================================
global VOX
VOX.IT = 1;                               % reset index

% get the onset of the next word
%----------------------------------------------------------------------
[Y,I,FS] = spm_voice_get_next(wfile);

% break if EOF
%----------------------------------------------------------------------
if isempty(I), error('Please repeat'); end


%% get 1 second segment and fundamental frequency
%======================================================================

% within Ockham's window W
%----------------------------------------------------------------------
LP   = log(P(:,1) + exp(-8));
W    = find(LP > (max(LP) - 3));

% get intervals (j) for this peak
%--------------------------------------------------------------------------
n    = numel(Y);
j    = fix((0:FS + FS/4) + I);
y    = Y(j(j < n & j > 1));
j    = logical(j < VOX.IT);
y(j) = 0;
J    = spm_voice_onsets(y,FS);

% get F0
%--------------------------------------------------------------------------
F0   = spm_voice_fundamental(Y,FS);

% likelihood search over fundamental and formant frequencies
%==========================================================================
F1    = exp(VOX.WHO(2).pE + VOX.R.F1);       % range of  formant frequencies
F     = zeros(numel(F1),1);
for i = 1:numel(F1)
    
    
    % set fundamental and formant frequencies
    %----------------------------------------------------------------------
    VOX.F0 = F0;
    VOX.F1 = F1(i);
    
    % sentence
    %----------------------------------------------------------------------
    clear xy
    for k = 1:1 % numel(J)
        xy(k,1) = spm_voice_ff(y(J{k}),FS);
    end
    
    % free energy
    %----------------------------------------------------------------------
    L    = spm_voice_likelihood(xy,W);
    L    = bsxfun(@plus,L,LP);
    Q    = spm_softmax(L);
    Fi   = sum(Q.*(L - log(Q + exp(-8))));
    F(i) = spm_softmax(Fi(:))'*Fi(:);
    
end

% posterior over F1
%--------------------------------------------------------------------------
F     = F - max(F);
L     = F(:);
L1    = F1(:);
L(:)  = spm_softmax(F);
F1    = L1'*spm_vec(sum(L,2));

% graphics if requested
%--------------------------------------------------------------------------
try VOX.formant;
    
    % plot
    %----------------------------------------------------------------------
    spm_figure('GetWin','Frequencies');
    subplot(3,1,3)
    plot(L1,F,'or',L1,F,':r',[F1 F1],[min(F) 0],'b')
    xlabel('frequency (Hz)'), ylabel('log likelihood')
    title('Formant frequency (F1)','FontSize',16), axis square, spm_axis tight
    drawnow
    
end


