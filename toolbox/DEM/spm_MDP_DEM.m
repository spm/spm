function DEM = spm_MDP_DEM(DEM,demi,O,o)
% auxiliary (link) function for mixed hierarchical (MDP/DEM) models
% FORMAT DEM = spm_MDP_DEM(DEM,demi,O,o)
%
% DEM      - DEM structure
% demi     - mapping from discrete outcomes to hidden causes
%  demi.C  - cell array of true causes for each combination of outcomes
%            the appropriate array is then placed in DEM.C
%  demi.U  - cell array of hidden causes for each combination of outcomes
%            the Bayesian model average is placed in DEM.U
% O{g}     - cell array of priors over discrete outcomes
% o(g x 1) - vector of true outcomes
%
% completes the following fields:
%   DEM.X{g} - posterior probability over g models and t times
%
% This routine performs a Bayesian model comparison using (DEM) Bayesian
% filtering and places the results in fields of the DEM structure; so that
% MDP schemes can pick them up as likelihood terms in the next hierarchical
% level. The outcomes of the (discrete) MDP scheme at the superordinate
% level specify the hidden causes at the current level. These enter as
% Bayesian model averages of the continuous causes. The resulting
% posterior over hidden causes then furnishes the posterior over outcomes
% using Bayesian model reduction, based on the free energy accumulated
% (integrated) over time. This free energy is supplemented with the prior
% over discrete outcomes; thereby constituting a posterior over outcomes.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_DEM.m 6901 2016-10-08 13:21:41Z karl $


% evaluate true and priors over causes given discrete states
%--------------------------------------------------------------------------
po    = spm_cross(O);
ind   = num2cell(o);
DEM.C = demi.C{ind{:}};
DEM.U = 0;
for i = 1:size(po,1)
    for j = 1:size(po,2)
        for k = 1:size(po,3)
            DEM.U = DEM.U + demi.U{i,j,k}*po(i,j,k);
        end
    end
end

% integrate system to generate data
%--------------------------------------------------------------------------
DEM.db = 0;
DEM    = spm_ADEM(DEM);

% posterior probability over discrete states models
%==========================================================================
nv    = size(DEM.U,1);
nt    = size(DEM.U,2);
ic    = (1:nv) + size(DEM.qU.C{1},1) - nv;
P     = spm_zeros(po);
for i = 1:ndims(P)
    x{i} = ones(size(P,i),1);
end

    
% prior potential
%--------------------------------------------------------------------------
F     = log(po);
    
% and accumulate log evidence
%--------------------------------------------------------------------------
for t = 1:nt
    
    % posteriors and priors for this time point
    %----------------------------------------------------------------------
    qE    = DEM.qU.v{end}(:,t);
    qC    = DEM.qU.C{t}(ic,ic);
    pE    = DEM.U(:,t);
    pC    = DEM.pU.C{t}(ic,ic);
    gt    = t > 1;
    for i = 1:size(po,1)
        for j = 1:size(po,2)
            for k = 1:size(po,3)
                rE       = demi.U{i,j,k}(:,end);
                F(i,j,k) = F(i,j,k) + gt*spm_log_evidence(qE,qC,pE,pC,rE,pC);
            end
        end
    end
    
    % Bayesian model comparison
    %----------------------------------------------------------------------
    P(:)  = spm_softmax(F(:));
    
    % marginal posteriors
    %----------------------------------------------------------------------
    for g = 1:numel(O)
        X{g}(:,t) = spm_dot(P,x,g);
    end
end

% return  probability over models (i.e., outcomes at subordinate level)
%--------------------------------------------------------------------------
DEM.X = X;







