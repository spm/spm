function [MDP] = spm_MDP_VB_sleep(MDP,BMR)
% Bayesian model reduction (sleep) for MDP models
% FORMAT [MDP] = spm_MDP_VB_sleep(MDP,BMR)
%
% MDP  - (inverted) MDP structure
%
% BMR.g - modality [default: 1]
% BMR.o - outcomes - that induce REM [default: {}]
% BMR.x - increase in concentration parameters for BMR [default: 8]
% BMR.f - hearing factors to sum over [default: 0]
% BMR.T - log Bayes factor threshold [default: 1/4]
% BMR.m - indicator function to enable BMR [@(i,i1,i2,i3,i4)1]
%
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
% (reduced) generative model.(c.f., REM sleep)
%
% See also: spm_MDP_log_evidence.m, spm_MDP_VB and spm_MDP_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_sleep.m 7679 2019-10-24 15:54:07Z spm $


% deal with a sequence of trials
%==========================================================================

% BMR options
%--------------------------------------------------------------------------
try, g   = BMR.g; catch, g = 1;   end
try, o   = BMR.o; catch, o = {};  end
try, x   = BMR.x; catch, x = 8;   end
try, f   = BMR.f; catch, f = 0;   end
try, T   = BMR.T; catch, T = 1/4; end

% model selection function
%--------------------------------------------------------------------------
if isfield(BMR,'m')
    m = BMR.m;
else
    m = @(i,i1,i2,i3,i4)1;
end

% Baysian model reduction - parameters
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    [sa,ra] = spm_MDP_VB_prune(MDP(end).a{g},MDP(1).a0{g},f,x,T,m);
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


function [sA,nA] = spm_MDP_VB_prune(qA,pA,f,x,T,m)
% FORMAT [sA,nA] = spm_MDP_VB_prune(qA,pA,f,x,T,m)
% qA - posterior expectations
% pA - prior expectations
% f  - hidden factor to integrate over [defult: 0]
% x  - prior counts [default: 8]
% T  - threshold for Bayesian model reduction [default: three]
%
% sA - reduced posterior expectations
% rA - reduced prior expectations
%__________________________________________________________________________

% defaults
%--------------------------------------------------------------------------
if nargin < 4, T = 3; end
if nargin < 3, x = 8; end

% identify noninformative subspaces
%--------------------------------------------------------------------------
ndim  = ndims(pA);
for i = 1:ndim
    d    = 1:ndim;
    d(i) = [];
    s    = pA < 4 & pA > 0;
    for j = 1:(ndim - 1)
        s = sum(s,d(j));
    end
    ind{i} = find(s(:));
end

% motifs of salient (low free free energy) abilities
%--------------------------------------------------------------------------
mot{1} = [1 0 0;0 1 0;0 0 1];
mot{2} = [1 0 0;0 0 1;0 1 0];
mot{3} = [0 1 0;1 0 0;0 0 1];
mot{4} = [0 1 0;0 0 1;1 0 0];
mot{5} = [0 0 1;1 0 0;0 1 0];
mot{6} = [0 0 1;0 1 0;1 0 0];

% model space: additional concentration parameters (i.e., precision)
%--------------------------------------------------------------------------
nmot  = numel(mot);
for i = 1:nmot*nmot;
    rA{i} = spm_zeros(pA);
end

for m1 = 1:numel(mot)
    [i1,i2] = find(mot{m1});
    for m2 = 1:numel(mot)
        i  = nmot*(m1 - 1) + m2;
        for j = 1:numel(i1)
            dA    = spm_zeros(pA);
            dA(1:3,i1(j),:,i2(j),4) = mot{m2};
            rA{i} = rA{i} + dA*x;
        end
    end
end

% score models using Bayesian model reduction
%--------------------------------------------------------------------------
for i = 1:numel(rA);
    G    = spm_MDP_log_evidence(qA,pA,pA + rA{i});
    F(i) = sum(G(isfinite(G)));
end

% find any model that has greater evidence than the parent model
%--------------------------------------------------------------------------
[F,j] = min(F);
if (F + T) < 0
    rA = rA{j};
else
    rA = spm_zeros(pA);
end

% reduced posterior and prior
%--------------------------------------------------------------------------
sA     = qA;
nA     = pA;
for i1 = 1:size(qA,2)
    for i2 = 1:size(qA,3)
        for i3 = 1:size(qA,4)
            for i4 = 1:size(qA,5)
                
                % get posteriors, priors and reduced priors
                %----------------------------------------------------------
                p  = pA(:,i1,i2,i3,i4);
                q  = qA(:,i1,i2,i3,i4);
                r  = rA(:,i1,i2,i3,i4);
                j  = find(p);
                p  = p(j);
                q  = q(j);
                r  = j(find(r(j)));
                
                % redistribute concentration parameters if indicated
                %----------------------------------------------------------
                if numel(r);
                    sA(:,i1,i2,i3,i4) = 0;
                    nA(:,i1,i2,i3,i4) = 0;
                    sA(r,i1,i2,i3,i4) = sum(q);
                    nA(r,i1,i2,i3,i4) = sum(p);
                else
                    sA(j,i1,i2,i3,i4) = q;
                    nA(j,i1,i2,i3,i4) = p;
                end
            end
        end
    end
end

return

% alternative formulation with specified indicator functions
%==========================================================================

% defaults
%--------------------------------------------------------------------------
if nargin < 5, m = @(i,i1,i2,i3,i4)1; end
if nargin < 4, T = 3; end
if nargin < 3, x = 8; end
if nargin < 2, f = 0; end

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
                F  = 0;
                if length(j) > 1
                    for i = 1:length(j);
                        if m(i,i1,i2,i3,i4)
                            r    = p;
                            r(i) = r(i) + x;
                            F(i) = spm_MDP_log_evidence(q,p,r);
                        else
                            F(i) = 16;
                        end
                    end
                end
                
                % eliminate parameter
                %----------------------------------------------------------
                [F,i] = min(F);
                mF(i1,i2,i3,i4) = F;
                iF(i1,i2,i3,i4) = j(i);
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
                if F < - T;
                    sA(:,i1,i2,i3,i4) = 0;
                    rA(:,i1,i2,i3,i4) = 0;
                    sA(i,i1,i2,i3,i4) = sum(q);
                    rA(i,i1,i2,i3,i4) = sum(p);
                else
                    sA(j,i1,i2,i3,i4) = q;
                    rA(j,i1,i2,i3,i4) = p;
                end
                
            end
        end
    end
end



