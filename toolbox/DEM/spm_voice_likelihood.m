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
% $Id: spm_voice_likelihood.m 7567 2019-04-04 10:41:15Z karl $

% defaults
%--------------------------------------------------------------------------
global VOX
if nargin < 2
    W = 1:size(VOX.LEX,1);
else
    W = W(:)';
end

% handle arrays
%==========================================================================
if numel(xY) > 1
    for i = 1:size(xY,1)
        for j = 1:size(xY,2)
            [Li,Mi,Ni] = spm_voice_likelihood(xY(i,j),W);
            L(:,:,i,j) = Li;
            M(:,:,i,j) = Mi;
            N(:,:,i,j) = Ni;
        end
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
ni    = numel(i);

% log likelihood over lexical outcomes
%==========================================================================
Q      = spm_vec(xY.Q) - spm_vec(VOX.Q);
L      = zeros(size(VOX.LEX)) - 1024;
method = 'likelihood';

switch method
    
    case {'likelihood'}
        
        % log likelihood - lexical
        %------------------------------------------------------------------
        for w = W
            for k = 1:size(VOX.LEX,2)
                E      = (Q(i) - VOX.LEX(w,k).qE(i));      % error
                L(w,k) = - E'*VOX.LEX(w,k).qP(i,i)*E/2;    % log likelihood
            end
        end
        
    case {'BMR'}
        
        % full posterior
        %------------------------------------------------------------------
        pE      = zeros(ni,1);
        pC      = speye(ni,ni)*16;
        P{1}.X  = speye(ni,ni);
        P{1}.C  = {speye(ni,ni)};
        P{2}.X  = pE;
        P{2}.C  = pC;
        C       = spm_PEB(Q(i),P,16,1);
        qE      = C{2}.E;
        qC      = C{2}.C;
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        for w = W
            for k = 1:size(VOX.LEX,2)
                rE     = VOX.LEX(w,k).qE(i);
                rC     = spm_inv(VOX.LEX(w,k).qP(i,i));
                L(w,k) = spm_log_evidence(qE,qC,pE,pC,rE,rC);
            end
        end
        
    case {'neuronal'}
        
        % gradient descent free energy
        %==================================================================
        ni    = 128;
        nw    = numel(VOX.LEX);
        for w = W
            E(:,w) = VOX.LEX(w).qP*(Q - VOX.LEX(w).qE);
        end
        
        
        % evidence accumulation
        %------------------------------------------------------------------
        v   = zeros(nw,1);
        s   = spm_softmax(v);
        tau = 1/4;
        for i = 1:ni
            qx = log(s);
            v  = E'*Q/2 - log(s);
            v  = v - mean(v);
            s  = spm_softmax(qx + v/tau);
            
            e(:,i) = v;
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
        
    end
end

return

