function [L,M,N] = spm_voice_likelihood(xY,LEX,PRO,WHO)
% returns the lexical likelihood
% FORMAT [L,M,N] = spm_voice_likelihood(xY,LEX,PRO,WHO)
%
% xY   - word    structure array
% LEX  - lexical structure array
% PRO  - prosody structure array
% WHO  - speaker structure array
%
% L    - log likelihood over lexicon
% M    - log likelihood over prodidy
% M    - log likelihood over speaker
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
% $Id: spm_voice_likelihood.m 7535 2019-03-03 20:27:09Z karl $


% handle arrays
%==========================================================================
if numel(xY) > 1
    for i = 1:numel(xY)
        [Li,Mi,Ni] = spm_voice_likelihood(xY(i),LEX,PRO,WHO);
        L(:,:,i)   = Li;
        M(:,:,i)   = Mi;
        N(:,:,i)   = Ni;
    end
    return
end

% defaults
%--------------------------------------------------------------------------
global voice_options
try, nu = voice_options.nu; catch, nu  = 8;   end  % order (formants)
try, nv = voice_options.nv; catch, nv  = 8;   end  % order (interval)

% precision wieghting
%--------------------------------------------------------------------------
Nu    = size(LEX(1).Q,1);
Nv    = size(LEX(1).Q,2);
nu    = min(nu,Nu);
nv    = min(nv,Nv);
wu    = sparse(1:nu,1,1,Nu,1);
wv    = sparse(1:nv,1,1,Nv,1);
W     = kron(wv,wu);

% log likelihood over lexical outcomes
%==========================================================================
Q     = spm_vec(xY.Q) - spm_vec(LEX(1).Q);
for w = 1:size(LEX,1)
    for k = 1:size(LEX,2)
        
        % log likelihood - lexical
        %------------------------------------------------------------------
        E      = (Q - LEX(w,k).qE).*W;            % error
        L(w,k) = - E'*(LEX(w,k).qP)*E/2;          % log likelihood
        
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
P     = spm_vec(xY.P) - spm_vec(PRO(1).P);
P     = P(PRO(1).i);
for p = 1:numel(PRO)
    for k = 1:numel(PRO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P(p) - PRO(p).pE(k);           % error
        M(k,p) = - E'*(PRO(p).pP(k))*E/2;       % log likelihood
        
        % add log prior ,if specified
        %------------------------------------------------------------------
        try
            M(k,p) = M(k,p) + PRO(p).pL(k);     % log posterior
        end
    end
end

if nargout < 3, return, end

% log likelihood over identity outcomes
%==========================================================================
P     = spm_vec(xY.P) - spm_vec(PRO(1).P);
P     = P(WHO(1).i);
for p = 1:numel(WHO)
    for k = 1:numel(WHO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P(p) - WHO(p).pE(k);           % error
        N(k,p) = - E'*(WHO(p).pP(k))*E/2;       % log likelihood
        
        % add log prior ,if specified
        %------------------------------------------------------------------
        try
            N(k,p) = N(k,p) + WHO(p).pL(k);     % log posterior
        end
    end
end

return

