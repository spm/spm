function DEM_demo_hierarchical_optmisation
% This is the same as spm_nlsi_GH but tries to model the free energy as a
% function of conditional expectations using a sparse mixture of scaled
% Gaussians. The objective is to account for local maxima when optimising
% free energy by recasting the problem in terms of a parameterised mapping 
% from conditional expectation to free energy explicitly.
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% set up model (a simple GLM)
%==========================================================================
ny   = 16;
np   = 2;
 
M.IS = inline('U*P','P','M','U');
M.pE = zeros(np,1);
M.pC = speye(np,np);
M.hE = 4;
 
% fixed parameters
%--------------------------------------------------------------------------
U    = randn(ny,np);
 
% free parameters
%--------------------------------------------------------------------------
P    = randn(np,1);
 
% data
%--------------------------------------------------------------------------
Y.y  = M.IS(P,M,U);
 
 
% invert (Laplace assumption)
%==========================================================================
[ep,cp,eh,f] = spm_nlsi_GN(M,U,Y);
 
% invert (hierarchical optimisation)
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN_H(M,U,Y);
