function [L,M,N] = spm_voice_likelihood(xY,w)
% Return the lexical likelihood
% FORMAT [L,M,N] = spm_voice_likelihood(xY,w)
%
% xY   - word structure array
% w    - indices of words in VOX.LEX to consider
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
%
% The likelihood can be estimated directly under the assumption of
% negligible random fluctuations on acoustic samples. Alternatively,
% parametric empirical Bayes (PEB) can be used to estimate observation
% noise, followed by Bayesian model reduction (BMR) to evaluate the
% (marginal) likelihood. In normal operation, the explicit likelihood
% scheme is used, with the opportunity to model the effects of (speech) in
% noise with an additional variable: VOX.noise (see main body of script).
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_voice_likelihood.m 7750 2019-12-05 17:54:29Z spm $


% defaults
%--------------------------------------------------------------------------
global VOX
nw    = numel(VOX.LEX);
if nargin < 2
    k = 1:nw;
else
    k = w(:)';
end

% handle arrays and cells
%==========================================================================

% many intervals, one contraint
%--------------------------------------------------------------------------
if numel(xY) > 1
    for i = 1:size(xY,1)
        for j = 1:size(xY,2)
            [Li,Mi,Ni] = spm_voice_likelihood(xY(i,j),k);
            L(:,i,j)   = Li;
            M(:,:,i,j) = Mi;
            N(:,:,i,j) = Ni;
        end
    end
    return
end

% many contraints, one interval
%--------------------------------------------------------------------------
if iscell(k)
    for i = 1:numel(k)
        [Li,Mi,Ni] = spm_voice_likelihood(xY,k{i});
        L(:,i)     = Li;
        M(:,:,i)   = Mi;
        N(:,:,i)   = Ni;
    end
    return
end


% precision wieghting, in terms of basis function coefficients to include
%--------------------------------------------------------------------------
jstr   = {'dur'};                                  % for lexical inference
pstr   = {'Tu','Tv','Tf','Tw'};                    % for prosody inference
j      = find(ismember({VOX.PRO.str},jstr));       % indices for prosody

% means and method
%--------------------------------------------------------------------------
W      = spm_vec(xY.W);
P      = spm_vec(xY.P);
R      = spm_vec(xY.R);
L      = zeros(numel(VOX.LEX),1)  - exp(16);
dP     = zeros(size(VOX.LEX(1).X,2) - 2,nw);
method = 'likelihood';

% log likelihood over lexical outcomes
%==========================================================================
switch method
    
    case {'likelihood'}
        
        % loop over words and variants
        %------------------------------------------------------------------
        for w = k
             
            % general linear model
            %--------------------------------------------------------------
            X    = VOX.LEX(w).X;                   % GLM
            qC   = VOX.LEX(w).qC;                  % covariance (lexical)
            pE   = VOX.LEX(w).pE;                  % expectation (pitch)
            pC   = VOX.LEX(w).pC;                  % covariance  (pitch)
            pP   = VOX.LEX(w).pP;                  % MAP prior           
            
            pP   = VOX.LEX(w).M*W + pP;            % parameters (pitch)
            Eq   = X*pP - W;                       % residuals (W)
            Ep   = pP - pE;                        % residuals (P)
            
            % save pitch parameters
            %--------------------------------------------------------------
            dP(:,w) = pP(3:end);
            
            % log posterior - lexical
            %--------------------------------------------------------------
            L(w) = - Eq'*(qC\Eq)/2 ...               % likelihood
                   - Ep'*(pC\Ep)/2 ...               % prior
                   + VOX.LEX(w).L;                   % entropy terms

            
            % shrink estimators in proportion to noise
            %--------------------------------------------------------------
            if isfield(VOX,'noise')
                L(w) = L(w)/VOX.noise;
            end
            
            % supplement with the likelihood of (duration) prodidy features
            %--------------------------------------------------------------
            E    = P(j) - VOX.LEX(w).dE(j);
            L(w) = L(w) - E'*(VOX.LEX(w).dC(j,j)\E)/2;
            
        end
        
        
    case {'BMR'}
        
        % full posterior
        %------------------------------------------------------------------
        ni      = numel(W);
        pE      = zeros(ni,1);
        pC      = speye(ni,ni)*16;
        B{1}.X  = speye(ni,ni);
        B{1}.C  = {speye(ni,ni)};
        B{2}.X  = pE;
        B{2}.C  = pC;
        C       = spm_PEB(W,B,16,1);
        qE      = C{2}.E;
        qC      = C{2}.C;
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        for w = k
            
            % log likelihood - lexical
            %--------------------------------------------------------------
            rE   = VOX.LEX(w).qE;
            rC   = VOX.LEX(w).qC;
            L(w) = spm_log_evidence(qE,qC,pE,pC,rE,rC);
            
            % supplement with the likelihood of duration
            %--------------------------------------------------------------
            E    = P(j) - VOX.LEX(w).dE(j);
            D    = -  E'*(VOX.LEX(w).dC(j,j)\E)/2;
            L(w) = L(w) + D;
            
        end
        
    otherwise
        disp('unknown method')
end


if nargout < 2, return, end


% update pitch parameters
%--------------------------------------------------------------------------
j     = find(ismember({VOX.PRO.str},pstr));
dP    = dP*spm_softmax(L);
P(j)  = P(j) + dP;

% log likelihood over prosody outcomes
%==========================================================================
P     = P - spm_vec(VOX.P);
for p = 1:numel(VOX.PRO)
    for k = 1:numel(VOX.PRO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = P(p) - VOX.PRO(p).pE(k);           % error
        M(k,p) = -  E'*(VOX.PRO(p).pC(k)\E)/2;      % log likelihood
        
    end
end

if nargout < 3, return, end

% update identity parameters
%--------------------------------------------------------------------------
R(1)  = R(1) - log(xY.P.inf(1));                    % ff0
R(2)  = R(2) + dP(3);                               % ff1

% log likelihood over identity outcomes
%==========================================================================
R     = R - spm_vec(VOX.R);
for p = 1:numel(VOX.WHO)
    for k = 1:numel(VOX.WHO(p).pE)
        
        % log likelihood
        %------------------------------------------------------------------
        E      = R(p) - VOX.WHO(p).pE(k);           % error
        N(k,p) = -  E'*(VOX.WHO(p).pC(k)\E)/2;      % log likelihood
        
    end
end
