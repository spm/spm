function [C] = spm_mc_loss_C(x,P)
% cost function for the mountain car problem
% problem
% FORMAT [C] = spm_mc_loss_C(x,P)
%
% x     - hidden states
% v     - exogenous inputs
% P.x,k - parameters for gradient function:     G(x(1),P.p)
% P.q,p - parameters for cost or loss-function: C(x(1),P.q)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mc_loss_C.m 3757 2010-03-08 11:41:53Z guillaume $
 
 
% gradient (G) (quadratic potential = (x(1) - P.x)^2*P.k/2)
%--------------------------------------------------------------------------
C   = abs(x(1,:) - 1) > 1/4;


