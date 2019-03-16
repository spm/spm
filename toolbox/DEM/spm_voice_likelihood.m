function [L,M,N] = spm_voice_likelihood(xY,W)
% returns the lexical likelihood
% FORMAT [L,M,N] = spm_voice_likelihood(xY,W)
%
% xY   - word    structure array
% W    - indices of words in VOX.LEX to consider
% 
% assumes the following structures are in the global structure VOX
% VOX.LEX  - lexical structure array
% VOX.PRO  - prosody structure array
% VOX.WHO  - speaker structure array
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
% $Id: spm_voice_likelihood.m 7545 2019-03-16 11:57:13Z karl $

% defaults
%--------------------------------------------------------------------------
global VOX
if nargin < 2, W = 1:size(VOX.LEX,1); end

% handle arrays
%==========================================================================
if numel(xY) > 1
    for i = 1:numel(xY)
        [Li,Mi,Ni] = spm_voice_likelihood(xY(i),W);
        L(:,:,i)   = Li;
        M(:,:,i)   = Mi;
        N(:,:,i)   = Ni;
    end
    return
end

% defaults
%--------------------------------------------------------------------------
try, nu = VOX.nu; catch, nu  = 8;   end  % order (formants)
try, nv = VOX.nv; catch, nv  = 8;   end  % order (interval)

% precision wieghting
%--------------------------------------------------------------------------
Nu    = size(VOX.Q,1);
Nv    = size(VOX.Q,2);
nu    = min(nu,Nu);
nv    = min(nv,Nv);
wu    = sparse(1:nu,1,1,Nu,1);
wv    = sparse(1:nv,1,1,Nv,1);
i     = find(kron(wv,wu));

% log likelihood over lexical outcomes
%==========================================================================
Q     = spm_vec(xY.Q) - spm_vec(VOX.Q);
for w = W
    for k = 1:size(VOX.LEX,2)
        
        % log likelihood - lexical
        %------------------------------------------------------------------
        E      = (Q(i) - VOX.LEX(w,k).qE(i));           % error
        L(w,k) = - E'*VOX.LEX(w,k).qP(i,i)*E/2;         % log likelihood
        
        % add log prior, if specified
        %------------------------------------------------------------------
        try
            L(w,k) = L(w,k) + VOX.LEX(w,k).pL(w,k);     % log posterior
        end
    end
end


if nargout < 2, return, end

% log likelihood over prosody outcomes
%==========================================================================
P     = spm_vec(xY.P) - spm_vec(VOX.P);
P     = P(VOX.i);
for p = 1:numel(VOX.PRO)
    for k = 1:numel(VOX.PRO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P(p) - VOX.PRO(p).pE(k);              % error
        M(k,p) = - E'*VOX.PRO(p).pP(k)*E/2;            % log likelihood
        
        % add log prior ,if specified
        %------------------------------------------------------------------
        try
            M(k,p) = M(k,p) + VOX.PRO(p).pL(k);        % log posterior
        end
    end
end

if nargout < 3, return, end

% log likelihood over identity outcomes
%==========================================================================
P     = spm_vec(xY.P) - spm_vec(VOX.P);
P     = P(VOX.j);
for p = 1:numel(VOX.WHO)
    for k = 1:numel(VOX.WHO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P(p) - VOX.WHO(p).pE(k);             % error
        N(k,p) = - E'*VOX.WHO(p).pP(k)*E/2;           % log likelihood
        
        % add log prior ,if specified
        %------------------------------------------------------------------
        try
            N(k,p) = N(k,p) + VOX.WHO(p).pL(k);       % log posterior
        end
    end
end

return

