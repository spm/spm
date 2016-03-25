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
% $Id: spm_MDP_VB_sleep.m 6756 2016-03-25 09:49:08Z karl $


% deal with a sequence of trials
%==========================================================================

% options
%--------------------------------------------------------------------------
try, g   = OPTIONS.g; catch, g = 1;  end
try, o   = OPTIONS.o; catch, o = {}; end
try, x   = OPTIONS.x; catch, x = 8;  end
try, f   = OPTIONS.f; catch, f = 0;  end

% Baysian model reduction - parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    [sa,ra] = spm_MDP_VB_prune(MDP(end).a{g},MDP(1).a0{g},f,x);
end


% reiterate expectation maximisation (or rapid eye movement sleep)
%--------------------------------------------------------------------------
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
    
    % Bayesian updating and updated parameters
    %----------------------------------------------------------------------
    REM    = spm_MDP_VB_X(REM);
    MDP.a  = REM(N).a;
    MDP.a0 = REM(N).a0;
    
else
    
    % otherwise, use reduced posteriors and priors
    %----------------------------------------------------------------------
    MDP.a{g}  = sa;
    MDP.a0{g} = ra;
end


function [sA,rA] = spm_MDP_VB_prune(qA,pA,f,x)
% FORMAT [sA,rA] = spm_MDP_VB_prune(qA,pA,f,x)
% qA - posterior expectations
% pA - prior expectations
% f  - hidden factor to integrate over [defult: 0]
% x  - prior counts [default: 8]
%
% sA - reduced posterior expectations
% rA - reduced prior expectations
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
if nargin < 4, x = 8; end
if nargin < 3, f = 0; end

% column-wise model comparison
%--------------------------------------------------------------------------
for i1 = 1:size(qA,2)
    for i2 = 1:size(qA,3)
        for i3 = 1:size(qA,4)
            for i4 = 1:size(qA,5)
                
                % get posteriors, priors and cycle over reduced priors
                %----------------------------------------------------------
                p  = pA(:,i1,i2,i3,i4);
                q  = qA(:,i1,i2,i3,i4);
                j  = find(p);
                p  = p(j);
                q  = q(j);
                
                % informative state?
                %----------------------------------------------------------
                for i = 1:length(j);
                    r    = p;
                    r(i) = x;
                    F(i) = spm_MDP_log_evidence(q,p,r);
                end
                
                % eliminate parameter
                %----------------------------------------------------------
                [F,i] = min(F);
                mF(i1,i2,i3,i4) = F;
                iF(i1,i2,i3,i4) = i;
            end
        end
    end
end

% pool over s{f}
%---------------------------------------------------------------------
if f, sF = sum(mF,f); end

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
                i  = iF(i1,i2,i3,i4);
                j  = find(p);
                p  = p(j);
                q  = q(j);
                
                % BMC
                %----------------------------------------------------------
                if f == 0
                    F = mF(i1,i2,i3,i4);
                elseif f == 1
                    F = sF( 1,i2,i3,i4);
                elseif f == 2
                    F = sF(i1, 1,i3,i4);
                elseif f == 3
                    F = sF(i1,i2, 1,i4);
                elseif f == 4
                    F = sF(i1,i2,i3, 1);
                end
                
                % eliminate parameter
                %----------------------------------------------------------
                if F < - 1;
                    sA(:,i1,i2,i3,i4) = 0;
                    rA(:,i1,i2,i3,i4) = 0;
                    sA(i,i1,i2,i3,i4) = sum(q(i));
                    rA(i,i1,i2,i3,i4) = 1;
                else
                    sA(j,i1,i2,i3,i4) = q;
                    rA(j,i1,i2,i3,i4) = p;
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


