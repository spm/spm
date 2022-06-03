function [MDP] = spm_MDP_VB_update(MDP,PDP,OPTIONS)
% Bayesian model reduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_update(MDP,PDP,OPTIONS)
%
% MDP  - MDP structure (prior to exposure)
% PDP  - MDP structure (after exposure)
%
% OPTIONS.d   - update of initial states of hidden factors d [default: []]
% OPTIONS.eta - forgetting rate [default: 1]
% OPTIONS.BMR - Bayesian model reduction options:
%
% BMR.g - Bayesian model reduction of modality g [default: 1]
% BMR.f - hidden factors to contract over [default: 0]
% BMR.o - outcomes - that induce REM [default: {}]
%
% MDP  - (updated) model structure: with updated MDP.a
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
% $Id: spm_MDP_VB_sleep.m 7679 2019-10-24 15:54:07Z spm $

% defaults
%--------------------------------------------------------------------------
try, OPTIONS.d;   catch, OPTIONS.d = [];  end
try, OPTIONS.eta; catch, OPTIONS.eta = 1; end

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
    try, BMR.g; catch, BMR.g = 1;  end       % outcome modality
    try, BMR.f; catch, BMR.f = 0;  end       % factors to contract over
    try, BMR.o; catch, BMR.o = {}; end       % outcomes for empirical BMS

    % Baysian model reduction - likelihood parameters
    %----------------------------------------------------------------------
    for g = BMR.g
        
        [qa,pa] = spm_MDP_VB_prune(PDP.a{g},MDP.a{g},MDP.C{g},BMR.f);
        
        
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
            REM       = spm_MDP_VB_XX(REM);
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

% forgetting parameter – eta
%--------------------------------------------------------------------------
eta = OPTIONS.eta;

% likelihood Dirichlet parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    for g = 1:numel(MDP.a)
        N = sum(MDP.a{g}(:));
        MDP.a{g} = PDP.a{g}*(N + eta*MDP.T)/(N + MDP.T);
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



% Bayesian model reduction subroutine
%==========================================================================
function [qA,pA] = spm_MDP_VB_prune(qA,pA,C,f)
% FORMAT [sA,rA] = spm_MDP_VB_prune(qA,pA,C,f)
% qA - posterior expectations
% pA - prior expectations
% C  - log preferences
% f  - hidden factor to contract over [default: 0]
%
% sA - reduced posterior expectations
% rA - reduced prior expectations
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
nd  = size(qA);                         % size of tensor

if nargin < 3, f = 0;              end  % no contraction
if nargin < 3, C = zeros(1,nd(1)); end  % no preferences

% some Dirichlet parameters over conditionally independent factors
%--------------------------------------------------------------------------
s   = prod(nd(f + 1));                  % contraction scaling
if f
    pA  = spm_sum(pA,f + 1);
    qA  = spm_sum(qA,f + 1);
end

% Bayesian model reduction (starting with the largest posterior)
%--------------------------------------------------------------------------
qA    = qA(:,:);
pA    = pA(:,:);

% gradients of expected information gain (i.e., expected free energy)
%--------------------------------------------------------------------------
rA    = pA;
n     = size(rA,2);
for i = 1:1
    
    % evaluate gradients of expected free energy
    %----------------------------------------------------------------------
    [E,dEdA] = spm_MDP_MI(rA,C);   
    dA       = rA + n*dEdA/16;
    j        = dA < 0 | ~isfinite(dA);
    dA(j)    = 0;
    
    % evaluate free energy (negative log evidence)
    %----------------------------------------------------------------------
    [F,sA]   = spm_MDP_log_evidence(qA,pA,dA);
    j        = F < -3;
    rA(:,j)  = dA(:,j);
    G(i)     = sum(F);
        
    % terminate if free energy converges
    %----------------------------------------------------------------------   
    if ~any(j)
        disp(i)
        break
    end
        
    % illustrate structure learning
    %----------------------------------------------------------------------
    if false
        
        spm_figure('GetWin','structure learning');
        %------------------------------------------------------------------
        subplot(2,2,1), imagesc(64*spm_norm(rA)), title('likelihood mapping'), axis image
        subplot(2,2,2), plot(rA), title('Prior Dirichlet counts'), axis square
        subplot(2,2,3), plot(G), title('reduced free energy'), axis square, drawnow
    end
    
end


% redistribute scaled parameters over contracted dimensions
%--------------------------------------------------------------------------
pA  = rA;
qA  = sA;
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