function [MDP] = spm_MDP_VB_sleep(MDP,OPTIONS)
% Bayesian model reduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_sleep(MDP,OPTIONS)
%
% MDP  - (inverted) MDP structure
%
% MDP  - (reduced) model structure: with reduced MDP.a
%
% This routine optimises the hyperparameters of a NDP model (i.e.,
% concentration parameters encoding probabilities. It uses Bayesian model
% reduction to evaluate the evidence for models with and without a
% particular parameter in the columns of MDP.a (c.f., SWS)
%
% If specified, the scheme will then recompute posterior beliefs about the
% model parameters based upon (fictive) outcomes generated under its
% (reduced) a generative model.(c.f., REM sleep)
%
% See also: spm_MDP_log_evidence.m, spm_MDP_VB and spm_MDP_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_VB_sleep.m 6748 2016-03-14 10:04:41Z karl $
 
 
% deal with a sequence of trials
%==========================================================================
 
% options
%--------------------------------------------------------------------------
try, g   = OPTIONS.g;   catch, g   = 1; end
try, REM = OPTIONS.REM; catch, REM = 0; end
 
% Baysian model reduction - parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    [sa,ra] = spm_MDP_VB_prune(MDP(end).a{g},MDP(1).a0{g});
end
 
 
% reiterate expectation maximisation (or rapid eye movement sleep)
%--------------------------------------------------------------------------
if REM
    
    % remove previous experience
    %----------------------------------------------------------------------
    N  = 64;
    A  = MDP.A;
    try, MDP = rmfield(MDP,'s'); end
    try, MDP = rmfield(MDP,'o'); end
    try, MDP = rmfield(MDP,'u'); end
    
    % and install a generative process and reset priors
    %----------------------------------------------------------------------
    % MDP.A{g}  = spm_norm(sa);
    MDP.a{g}  = ra;
    MDP.a0{g} = ra;
    
    for i = 1:N
        MDP(i) = MDP(1);
    end
    
    % Bayesian updating and updated parameters
    %----------------------------------------------------------------------
    MDP   = spm_MDP_VB_X(MDP);
    MDP   = MDP(N);
    MDP.A = A;
    
else
    
    % otherwise, use reduced posteriors and priors
    %----------------------------------------------------------------------
    MDP.a{g}  = sa;
    MDP.a0{g} = ra;
end
 
 
function [sA,rA] = spm_MDP_VB_prune(qA,pA)
% FORMAT [sA,rA] = spm_MDP_VB_prune(qA,pA)
% qA - posterior expectations
% pA - prior expectations
% sA - reduced posterior expectations
% rA - reduced prior expectations
%__________________________________________________________________________

% column-wise reduction
%--------------------------------------------------------------------------
sA     = qA;
rA     = pA;
for i1 = 1:size(qA,2)
    for i2 = 1:size(qA,3)
        for i3 = 1:size(qA,4)
            for i4 = 1:size(qA,5)
                
                % get posteriors, priors and cycle over reduced priors
                %----------------------------------------------------------
                p  = pA(:,i1,i2,i3,i4);
                q  = qA(:,i1,i2,i3,i4);
                for i = 1:size(qA,1);
                    if p(i)
                        r    = p;
                        r(i) = 1/2;
                        dF   = spm_MDP_log_evidence(q,p,r);
                        
                        % eliminate parameter
                        %--------------------------------------------------
                        if dF < - 1/128;
                            sA(i,i1,i2,i3,i4) = 0;
                            rA(i,i1,i2,i3,i4) = 0;
                        end
                    else
                        sA(i,i1,i2,i3,i4) = 0;
                        rA(i,i1,i2,i3,i4) = 0;
                    end
                end
            end
        end
    end
end
 
 
function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                s = sum(A(:,i,j,k,l),1);
                if s, A(:,i,j,k,l) = A(:,i,j,k,l)/s; end
            end
        end
    end
end


