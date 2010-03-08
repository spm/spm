function [f] = spm_mc_fx(x,v,P)
% equations of motion for the mountain car problem using basis functions
% problem
% FORMAT [f] = spm_mc_fx(x,v,P)
%
% x     - hidden states
% v     - exogenous inputs
% P.x,k - parameters for gradient function:     G(x(1),P.p)
% P.q   - parameters for cost or loss-function: C(x(1),P.q)
%
% returns f = dx/dt = f  = [x(2);
%                           G - x(2)*C(x(1))]*dt;
%
% where C determines divergence of flow x(2) at any position x(1).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mc_fx.m 3757 2010-03-08 11:41:53Z guillaume $
 
 
% gradient (G)
%--------------------------------------------------------------------------
G   = spm_mc_loss_G(x(1),P);
 
% cost function (C)
%--------------------------------------------------------------------------
c   = spm_mc_loss_C(x(1),P);
 
% flow
%--------------------------------------------------------------------------
dt  = 1/4;
C   = P.p - P.q*(1 - c);
f   = [x(2); -G*c + x(2)*C]*dt;


