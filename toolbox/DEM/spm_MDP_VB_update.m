function [MDP] = spm_MDP_VB_update(MDP,PDP,OPTIONS)
% Bayesian model reduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_update(MDP,PDP,OPTIONS)
%
% MDP  - MDP structure (prior to exposure)
% PDP  - MDP structure (after exposure)
%
% OPTIONS.d   - update of initial states of hidden factors d [default: []]
% OPTIONS.BMR - Bayesian model reduction options:
%
% BMR.g - Bayesian model reduction of modality g [default: all]
% BMR.f - hidden factors to contract over [default: 0]
% BMR.o - outcomes - that induce REM [default: {}]
% BMR.T - Occams threshold [default: 2]
%
% MDP  - (updated) model structure: with updated MDP.a
%
% This routine optimises the hyperparameters of a POMDP model (i.e.,
% concentration parameters encoding likelihoods). It uses Bayesian model
% reduction to evaluate the evidence for models with and without an changes
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
% See also: spm_MDP_log_evidence.m, spm_MDP_VB and spm_MDP_VB_sleep.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% defaults
%--------------------------------------------------------------------------
try, OPTIONS.d;   catch, OPTIONS.d   = []; end

% deal with multiple agents
%==========================================================================
if numel(MDP) > 1
    for i = 1:numel(MDP)
        MDP(i) = spm_MDP_VB_update(MDP(i),PDP(i),OPTIONS);
    end
    return
end

% update initial states (post-diction)
%==========================================================================
for f = OPTIONS.d
    MDP.D{f} = PDP.X{f}(:,end);
end

% Bayesian model reduction
%==========================================================================
if isfield(OPTIONS,'BMR')
    
    % BMR options
    %----------------------------------------------------------------------
    BMR = OPTIONS.BMR;
    try, BMR.g; catch, BMR.g = 0;  end       % outcome modality
    try, BMR.f; catch, BMR.f = 0;  end       % factors to contract over
    try, BMR.o; catch, BMR.o = {}; end       % outcomes for empirical BMS
    try, BMR.T; catch, BMR.T = 0;  end       % outcomes for empirical BMS

    % default: all modalities
    %----------------------------------------------------------------------
    if ~BMR.g, BMR.g = 1:numel(MDP.a); end

    % Baysian model reduction - likelihood parameters
    %----------------------------------------------------------------------
    for g = BMR.g
        
        % structure learning with BMR
        %==================================================================
        [qa,pa] = spm_MDP_VB_prune(PDP.a{g},MDP.a{g},BMR.f,BMR.T,MDP.C{g});
                
        % Dirichlet accumulation under reduced model and outcomes
        % (c.f., consolidation during rapid eye movement sleep)
        %==================================================================
        N  = numel(BMR.o);
        if N
            
            % remove previous experience
            %--------------------------------------------------------------
            REM      = MDP;
            try, REM = rmfield(REM,'s'); end
            try, REM = rmfield(REM,'o'); end
            try, REM = rmfield(REM,'u'); end
            
            % and install a generative process and reset priors
            %--------------------------------------------------------------
            REM.a{g}  = pa;
            REM.o = BMR.o{1};
            for i = 1:N
                REM(1,i)   = REM(1);
                REM(1,i).o = BMR.o{i};
            end
            
            % Bayesian belief updating: posteriors and priors
            %--------------------------------------------------------------
            REM       = spm_MDP_VB_XXX(REM);
            PDP.a     = REM(N).a;
            
        else
            
            % otherwise, use reduced posteriors and priors
            %--------------------------------------------------------------
            PDP.a{g}  = qa;
            
        end
    end
end

% update Dirichlet parameters: move Dirichlet parameters from PDP to MDP
%==========================================================================

% likelihood Dirichlet parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    for g = 1:numel(MDP.a)
        MDP.a{g} = PDP.a{g};
    end
end

% check for remaining concentration parameters at this level
%--------------------------------------------------------------------------
try,  MDP.b = PDP.b; end
try,  MDP.c = PDP.c; end
try,  MDP.d = PDP.d; end
try,  MDP.e = PDP.e; end

% check for concentration parameters at nested levels
%--------------------------------------------------------------------------
try,  MDP.MDP(1).a = PDP.mdp(end).a; end
try,  MDP.MDP(1).b = PDP.mdp(end).b; end
try,  MDP.MDP(1).c = PDP.mdp(end).c; end
try,  MDP.MDP(1).d = PDP.mdp(end).d; end
try,  MDP.MDP(1).e = PDP.mdp(end).e; end

return
