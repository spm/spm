function [DEM] = spm_MDP_DEM(DEM,pP,p)
% auxiliary (link) function for mixed hierarchical (MDP/DEM) models
% FORMAT DEM = spm_MDP_DEM(DEM,pP,p)
%
% DEM - DEM structure
% pP  - prior probabilities over model parameters
% pP  - index of true parameters set for generative process
%
% expects the following fields:
%   DEM.P - a structure array of alternative parameterisations
%
% completes the following fields:
%   DEM.X - posterior probability over models (i.e., parameterisations)
%   DEM.s - true parameterisation
%
% This routine performs a Bayesian model comparison using (DCM) Bayesian
% filtering and places the results in fields of the Dem structure; so that
% MDP schemes can pick them up as likelihood terms in the next hierarchical
% models. The generative process for the DCM model is assumed to be the
% same – and models aare defined in terms of their parameterisation.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_DEM.m 6866 2016-09-05 09:19:42Z karl $


% integrate system to generate data
%--------------------------------------------------------------------------
DEM.M.pE = DEM.P{p};
DEM.G.pE = DEM.P{p};
DEM      = spm_ADEM(DEM);

% revisit outcomes under alternative models
%--------------------------------------------------------------------------
for i = 1:numel(pP)
    DEM.M.pE = DEM.P{i};
    dem(i)   = spm_DEM(DEM);
    F(i,1)     = dem(i).F; 
end

% Bayesian model comparison
%--------------------------------------------------------------------------
P     = spm_softmax(F);

% return best model and probability over models
%--------------------------------------------------------------------------
[i,j] = max(F);
DEM   = dem(j);
DEM.X = P;
DEM.s = p;





