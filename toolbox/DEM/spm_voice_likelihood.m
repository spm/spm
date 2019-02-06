function [L,M] = spm_voice_likelihood(xY,LEX,PRO)
% returns the lexical likelihood
% FORMAT [L,M] = spm_voice_likelihood(xY,LEX,PRO)
%
% xY   - word    structure array
% LEX  - lexical structure array
% PRO  - prosody structure array
%
% L    - log likelihood over lexicon
% M    - log likelihood over prodidy
%
% This routine returns the log likelihood of a word and prosody based upon
% a Gaussian mixture model; specified in terms of a prior expectation and
% precision for each word (or prosody).  Prosody is categorised over
% several  dimensions (i.e., eigenmodes).  For both lexical and prosody,
% likelihoods are evaluated based upon the deviations from the expected
% parameters, over all words and prosody dimensions.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_likelihood.m 7528 2019-02-06 19:19:49Z karl $



% handle arrays
%==========================================================================
if numel(xY) > 1
    for i = 1:numel(xY)
        [Li,Mi]  = spm_voice_likelihood(xY(i),LEX,PRO);
        L(:,:,i) = Li;
        M(:,:,i) = Mi;
    end
    return
end

% log likelihood over lexical outcomes
%==========================================================================
Q     = spm_vec(xY.Q)- spm_vec(LEX(1).Q);
Q     = LEX(1).U'*Q;
for w = 1:size(LEX,1)
    for k = 1:size(LEX,2)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = Q - LEX(w,k).pE;               % error
        L(w,k) = - E'*(LEX(w,k).pP)*E/2;        % log likelihood

        % add log prior, if specified
        %------------------------------------------------------------------
        try
            L(w,k) = L(w,k) + LEX(w,k).pL(w,k); % log posterior
        end
    end
end


if nargout < 2, return, end

% log likelihood over prosody outcomes
%==========================================================================
for p = 1:numel(PRO)
    P = spm_vec(xY.P) - spm_vec(PRO(1).P);
    P = PRO(1).U'*P;
    for k = 1:numel(PRO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P - PRO(p).pE(k);              % error
        M(k,p) = - E'*(PRO(p).pP(k))*E/2;       % log likelihood
        
        % add log prior ,if specified
        %------------------------------------------------------------------
        try
            M(k,p) = M(k,p) + PRO(p).pL(k);     % log posterior
        end
    end
end

return

% add phonetic priors on lexical liklihood
%==========================================================================
if isfield(PRO,'J')
    J    = PRO(1).J;
    LP   = spm_dot(J,num2cell(spm_softmax(M),1));
    L    = bsxfun(@plus,L,log(LP + exp(-4)));
end
