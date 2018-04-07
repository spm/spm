function DCM = DEMO_AI_NLSI
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_AI_NLSI.m 7287 2018-04-07 13:15:09Z karl $
 
% set up and preliminaries
%==========================================================================
% rng('default')
  
% specify generative model
%--------------------------------------------------------------------------
n     = 32;
m     = 1;
M.X   = randn(n,m);
M.pE  = zeros(m,1);
M.pC  = speye(m,m)*32;
M.hE  = zeros(1,1);
M.hC  = speye(1,1);
M.Q   = speye(n,n);
P.p   = 2;%randn(m,1);
P.h   = -1;%randn(1,1);

% generate data
%--------------------------------------------------------------------------
Y     = M.X*P.p + spm_sqrtm(spm_inv(exp(P.h)*M.Q))*randn(n,1);

% log likelihood (energy) function ln P(Y,P)
%--------------------------------------------------------------------------
M.L   = @(P,Y,M)spm_logdet(exp(P.h)*M.Q)/2 - ...
        exp(P.h)*(M.X*P.p - Y)'*M.Q*(M.X*P.p - Y)/2 - ...
        (M.pE - P.p)'*(M.pC\(M.pE - P.p))/2;
    
% invert
%--------------------------------------------------------------------------
spm_nlsi_AI(M,Y);