% Bayesian model reduction subroutine
%==========================================================================
function [qA,pA] = spm_MDP_VB_prune(qA,pA,f,T,C,OPT)
% FORMAT [sA,rA] = spm_MDP_VB_prune(qA,pA,f,T,C,OPT)
% qA  - posterior expectations
% pA  - prior expectations
% T   - threshold for Occam's window (natural units) [default: 2]
% f   - hidden factor to contract over [default: 0]
% C   - log preferences
% OPT - 'MI' or 'SIMPLE'
%
% sA  - reduced posterior expectations
% rA  - reduced prior expectations
%
% OPT - 'MI': mutual information
%
% This implementation compares models (i.e., prior concentration
% parameters) that change in the direction of maximising expected free
% energy; namely, maximising the mutual information entailed by a
% likelihood mapping or transition matrix (plus the log preference over
% outcomes or states). If the reduced prior exceeds the specified Occams
% window, in terms of the reduced free energy, the reduced  priors and
% posteriors replace the full priors and posteriors. Effectively, this
% implements the structural hyperprior that likelihood mappings with a high
% mutual information are plausible and accepts these new priors if there is
% sufficient evidence for them. This can be regarded as a generic form of
% structure learning.

% OPT - 'SIMPLE':  simple suppression of implausible likelihoods
%
% This option suppresses entries of the likelihood tensor with small values
% and places the posterior probability mass over more likely outcomes. In
% effect,  it endows the first dimension of the likelihood tensor with a
% sparsity by setting elements with low Dirichlet counts to 0; thereby
% precluding any further evidence accumulation.
%
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
nd  = size(qA);                           % size of tensor
if nargin < 3, f   = 0;              end  % no contraction
if nargin < 4, T   = 0;              end  % Occams threshold
if nargin < 5, C   = zeros(1,nd(1)); end  % no preferences
if nargin < 6, OPT = 'MI';           end  % BMR type

s   = prod(nd(f + 1));                    % contraction scaling

% sum Dirichlet parameters over conditionally independent factors
%--------------------------------------------------------------------------
if f
    pA = spm_sum(pA,f + 1);
    qA = spm_sum(qA,f + 1);
end

% unfold tensors
%--------------------------------------------------------------------------
qA  = qA(:,:);
pA  = pA(:,:);

% gradients of expected information gain (i.e., expected free energy)
%==========================================================================

switch OPT

    case('MI')

        % evaluate gradients of expected free energy
        %------------------------------------------------------------------
        [E,dEdA] = spm_MDP_MI(qA,C);
        rA       = qA.*exp(512*qA.*dEdA);       % reduced prior
        rA       = times(rA,sum(qA)./sum(rA));  % preserve probability mass


    case('SIMPLE')

        % retain reduced priors and posteriors if outwith Occam's window
        %------------------------------------------------------------------
        rA       = spm_psi(qA);                 % reduced (log) prior
        rA       = minus(rA,max(rA));           % normalise
        rA       = times(qA,rA > -2);           % prune
        rA       = times(rA,sum(qA)./sum(rA));  % preserve probability mass


    otherwise

end

% retain reduced priors and posteriors if outwith Occam's window
%--------------------------------------------------------------------------
F        = spm_MDP_log_evidence(qA,pA,rA);
j        = F < T;
qA       = pA;
qA(:,j)  = rA(:,j);
pA(:,j)  = rA(:,j);


% redistribute scaled parameters over contracted dimensions
%--------------------------------------------------------------------------
if f
    pA = plus(zeros(nd),pA/s);
    qA = plus(zeros(nd),qA/s);

else

    % fold tensors
    %----------------------------------------------------------------------
    qA       = reshape(qA,nd);
    pA       = reshape(pA,nd);

end



return
