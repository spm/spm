function [G] = spm_mc_loss_G(x,P)
% assumed potential gradient for the mountain car problem
% problem
% FORMAT [G] = spm_mc_loss_G(x,P)
%
% x     - hidden states
% v     - exogenous inputs
% P.x,k - parameters for gradient function:     G(x(1),P.p)
% P.q,q - parameters for cost or loss-function: C(x(1),P.q)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
 
% gradient (G) (quadratic potential = (x(1) - P.x)^2*P.k/2)
%--------------------------------------------------------------------------
G   = (x(1,:) - P.x)*P.k;


