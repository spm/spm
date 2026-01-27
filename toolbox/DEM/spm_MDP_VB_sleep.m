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
% BMR.T - log Bayes factor threshold [default: 0]
%
% MDP  - (reduced) model structure: with reduced MDP.a
%
% This routine optimises the hyperparameters of a POMDP model (i.e.,
% concentration parameters encoding likelihoods). It uses Bayesian model
% reduction to evaluate the evidence for models with and without changes
% in Dirichlet counts (c.f., SWS or introspection)
%
% If specified, the scheme will then recompute posterior beliefs about the
% model parameters based upon (fictive) outcomes generated under its
% (reduced) generative model.(c.f., REM sleep)
%
% This version compares models (i.e., prior concentration parameters) that
% change in the direction of maximising expected free energy; namely,
% maximising the mutual information entailed by a likelihood mapping or
% transition matrix (plus the log preference over outcomes or states). If
% the reduced prior exceeds the specified Occams window, in terms of the
% reduced free energy, the reduced  priors and posteriors replace the full
% priors and posteriors. Effectively, this implements the structural
% hyperprior that likelihood mappings with a high mutual information are
% plausible and accepts these new priors if there is sufficient evidence
% for them. This can be regarded as a generic form of structure learning.
%
% See also: spm_MDP_log_evidence.m, spm_MDP_VB and spm_MDP_VB_update.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% deal with a sequence of trials
%==========================================================================
if ~isfield(MDP,'a'), return, end

% BMR options
%--------------------------------------------------------------------------
try, g  = BMR.g; catch, g = 1;    end       % outcome modality
try, o  = BMR.o; catch, o = {};   end       % outcomes for empirical BMS
try, f  = BMR.f; catch, f = 0;    end       % factors to contract over
try, T  = BMR.T; catch, T = 0;    end       % threshold or window (nats)

% for multiple agents
%--------------------------------------------------------------------------
if size(MDP,1) > 1
    for i = 1:size(MDP,1)
        MDP(i) = spm_MDP_VB_sleep(MDP(i),BMR);
    end
    return
end

% structure learning with BMR
%==========================================================================
for i = 1:numel(g)

    % Prior Dirichlet counts [default: uninformative and unitary]
    %----------------------------------------------------------------------
    if isfield(MDP,'a0')
        a0 = MDP.a0{g(i)};
    else
        a0 = (MDP.a{g(i)} > 0);
    end

    [s,r] = spm_MDP_VB_prune(MDP.a{g(i)},a0,f,T);
    sa{i} = s;
    ra{i} = r;

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
    for i = 1:numel(g)
        REM.a{g(i)}  = ra{i};
        REM.a0{g(i)} = ra{i};
    end

    REM.o = o{1};
    for i = 1:N
        REM(i)   = REM(1);
        REM(i).o = o{i};
    end
    
    % Bayesian belief updating: posteriors and priors
    %----------------------------------------------------------------------
    REM       = spm_MDP_VB_XXX(REM);
    MDP.a     = REM(N).a;
    MDP.a0{g} = REM(N).a0;
    
else
    
    % otherwise, use reduced posteriors and priors
    %----------------------------------------------------------------------
    for i = 1:numel(g)
        MDP.a{g(i)}  = sa{i};
        MDP.a0{g(i)} = ra{i};
    end

end

return
