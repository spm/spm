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
% $Id: spm_voice_identity.m 7597 2019-05-23 18:42:38Z karl $


%% log prior over lexical content
%==========================================================================
global VOX
VOX.IT = 1;                               % reset index

% get F0
%----------------------------------------------------------------------
[Y,I,FS] = spm_voice_get_next(wfile);
F0       = spm_voice_fundamental(Y,FS);
rE       = 28 + F0/20;
rC       = 4;
dR       = linspace(-2,2,8)*sqrt(rC);
R1       = dR + rE;
LR       = -dR.^2/rC/2;


%% deal with multiple words
%--------------------------------------------------------------------------
for w = 1:size(P,2)
    
    % get the onset of the next word
    %----------------------------------------------------------------------
    [Y,I,FS] = spm_voice_get_next(wfile);
    
    % break if EOF
    %----------------------------------------------------------------------
    if isempty(I), 
        if w == 1
            error('Please repeat'); 
        else
            break
        end
    end
    
    % get 1 second segment and fundamental frequency
    %======================================================================
    
    % within Ockham's window W
    %----------------------------------------------------------------------
    LP    = log(P(:,w) + exp(-8));
    k     = find(LP > (max(LP) - 3));
    
    % get intervals (j) for this peak
    %----------------------------------------------------------------------
    n     = numel(Y);
    j     = fix((0:FS + FS/4) + I);
    y     = Y(j(j < n & j > 1));
    j     = logical(j < VOX.IT);
    y(j)  = 0;
    J     = spm_voice_onsets(y,FS);
    
    % likelihood search over fundamental and formant frequencies
    %======================================================================
    LL    = zeros(numel(R1),1);
    for i = 1:numel(R1)
        
        % set fundamental and formant frequencies
        %------------------------------------------------------------------
        VOX.F0 = F0;
        VOX.F1 = R1(i);
        
        % sentence
        %------------------------------------------------------------------
        clear xy
        for k = 1:numel(J)
            xy(k,1) = spm_voice_ff(y(J{k}),FS);
        end
        
        % free energy
        %------------------------------------------------------------------
        L     = spm_voice_likelihood(xy,k);
        L     = bsxfun(@plus,L,LP);
        Q     = spm_softmax(L);
        F     = sum(Q.*(L - log(Q + exp(-8))));
        [F,m] = max(F);
        LL(i) = F;
        
    end
    
    % advance pointer and accumulate evidence
    %----------------------------------------------------------------------
    VOX.IT  = I + fix(J{m}(end));
    FF(:,w) = LL;
    
end

% posterior over F1
%--------------------------------------------------------------------------
P     = spm_softmax(sum(FF,2) + LR(:));
F     = log(P + exp(-64));
F1    = R1*P;

% graphics if requested
%--------------------------------------------------------------------------
try VOX.formant;
    
    % plot
    %----------------------------------------------------------------------
    spm_figure('GetWin','Frequencies');
    subplot(3,2,6)
    plot(R1,F,'or',R1,F,':r',[F1 F1],[min(F) 0],'b')
    xlabel('frequency (Hz)'), ylabel('log likelihood')
    title('Formant frequency (F1)','FontSize',16)
    axis square, spm_axis tight
    drawnow
    
end


