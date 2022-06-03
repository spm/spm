function [MDP] = spm_MDP_VB_sleep(MDP,BMR)
% Bayesian model reduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_sleep(MDP,BMR)
%
% MDP  - (inverted) MDP structure
%
% BMR.g - modality [default: 1]
% BMR.o - outcomes - that induce REM [default: {}]
% BMR.x - increase in concentration parameters for BMR [default: 8]
% BMR.f - hidden factors to contract over [default: 0]
% BMR.T - log Bayes factor threshold [default: 2]
%
% MDP  - (reduced) model structure: with reduced MDP.a
%
% This routine optimises the hyperparameters of a MDP model (i.e.,
% concentration parameters encoding likelihoods. It uses Bayesian model
% reduction to evaluate the evidence for models with and without an
% increase in a particular parameter in the columns of MDP.a (c.f., SWS)
%
% If specified, the scheme will then recompute posterior beliefs about the
% model parameters based upon (fictive) outcomes generated under its
% (reduced) generative model.(c.f., REM sleep)
%
% This version compares models (i.e., prior concentration parameters) in
% which one unique outcome is much more likely than any other outcome.
% Furthermore, it will then reduce posterior concentration parameters to
% prior concentration parameters for all likelihood mappings, to and from
% the parameter in question, if the reduced prior exceeds the specified
% Occams window. effectively, this implements the structural hyperprior
% that likelihood mappings with a high mutual information are plausible.
%
% See also: spm_MDP_log_evidence.m, spm_MDP_VB and spm_MDP_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_sleep.m 8262 2022-06-03 14:15:28Z karl $


% deal with a sequence of trials
%==========================================================================

% BMR options
%--------------------------------------------------------------------------
try, g  = BMR.g; catch, g = 1;    end       % outcome modality
try, o  = BMR.o; catch, o = {};   end       % outcomes for empirical BMS
try, f  = BMR.f; catch, f = 0;    end       % factors to contract over
try, T  = BMR.T; catch, T = 2;    end       % threshold or window (nats)

% for multiple agents
%--------------------------------------------------------------------------
if size(MDP,1) > 1
    for i = 1:size(MDP,1)
        MDP(i) = spm_MDP_VB_sleep(MDP(i),BMR);
    end
    return
end

% Prior Dirichlet counts (assumed to be one, everywhere)
%--------------------------------------------------------------------------
if isfield(MDP,'a0')
    a0 = MDP.a0{g};
else
    a0 = (MDP.a{g} > 0);
end

% Baysian model reduction - likelihood parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    [sa,ra] = spm_MDP_VB_prune(MDP.a{g},a0,f,T);
else
    return
end

% repeat active inference and learning under reduced model and outcomes
% (c.f., rapid eye movement sleep)
%==========================================================================
N  = numel(o);
if N
    
    % remove previous experience
    %----------------------------------------------------------------------
    REM  = MDP;
    try, REM = rmfield(REM,'s'); end
    try, REM = rmfield(REM,'o'); end
    try, REM = rmfield(REM,'u'); end
    
    % and install a generative process and reset priors
    %----------------------------------------------------------------------
    REM.a{g}  = ra;
    REM.a0{g} = ra;
    REM.o = o{1};
    for i = 1:N
        REM(i)   = REM(1);
        REM(i).o = o{i};
    end
    
    % Bayesian belief updating: posteriors and priors
    %----------------------------------------------------------------------
    REM       = spm_MDP_VB_X(REM);
    MDP.a     = REM(N).a;
    MDP.a0{g} = ra;
    
else
    
    % otherwise, use reduced posteriors and priors
    %----------------------------------------------------------------------
    MDP.a{g}  = sa;
    MDP.a0{g} = ra;
    
end



return



% Bayesian model reduction subroutine
%==========================================================================
function [qA,pA] = spm_MDP_VB_prune(qA,pA,f,T)
% FORMAT [sA,rA] = spm_MDP_VB_prune(qA,pA,f,T)
% qA - posterior expectations
% pA - prior expectations
% f  - hidden factor to contract over [default: 0]
% T  - Occam's window [default: 2]
%
% sA - reduced posterior expectations
% rA - reduced prior expectations
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
if nargin < 3, f = 0; end
if nargin < 4, T = 2; end

% some Dirichlet parameters over conditionally independent factors
%--------------------------------------------------------------------------
nd  = size(pA);                         % size of likelihood tensor
s   = prod(nd(f + 1));                  % contraction scaling
if f
    pA  = spm_sum(pA,f + 1);
    qA  = spm_sum(qA,f + 1);
end

% Bayesian model reduction (starting with the largest posterior)
%--------------------------------------------------------------------------
qA    = qA(:,:);
pA    = pA(:,:);
c     = ones(1,size(pA,2));
while(any(c))
    
    % find maximum Dirichlet parameter for this hidden state
    %----------------------------------------------------------------------
    [d,j]  = max(max(qA - pA).*c);
    [d,i]  = max(qA(:,j));
    c(j)   = 0;
    
    % create reduced prior (by adding eight Dirichlet counts)
    %----------------------------------------------------------------------
    rA     = pA(:,j);
    rA(i)  = rA(i)*4*s;

    % evaluate reduced posterior and free energy (negative log evidence)
    %----------------------------------------------------------------------
    [F,sA] = spm_MDP_log_evidence(qA(:,j),pA(:,j),rA);
    
    % retain reduced priors and posteriors if outwith Occam's window
    %----------------------------------------------------------------------
    if F < -T
        disp('mutual'),disp(sum(F))
        qA(i,:) = pA(i,:);            % forget afferent mapping (column)
        qA(:,j) = pA(:,j);            % forget efferent mapping (row)
        qA(:,j) = sA;
        pA(:,j) = rA;
    end
end

% redistribute scaled parameters over contracted dimensions
%--------------------------------------------------------------------------
if f
    pA = bsxfun(@plus,zeros(nd),pA/s);
    qA = bsxfun(@plus,zeros(nd),qA/s);
end

return


function A  = spm_log(A)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
A           = log(A);
A(isinf(A)) = -32;

function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);