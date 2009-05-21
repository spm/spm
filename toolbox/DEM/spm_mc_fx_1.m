function [f] = spm_mc_fx_1(x,v,P)
% equations of motion for the mountain car problem using basis functions
% problem
% FORMAT [f] = spm_mc_fx_1(x,v,P)
%
% x   - hidden states
% v   - exogenous inputs
% P.p - parameters for gradient function:     G(x(1),P.p)
% P.q - parameters for cost or loss function: C(x(1),P.q)
%
% returns f = dx/dt = f  = [x(2);
%                           G - x(2)*C]*dt;
%
% where C determines divergence of flow x(2) at any position x(1).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mc_fx_1.m 3140 2009-05-21 18:38:17Z karl $
 
 
% gradient (G)
%--------------------------------------------------------------------------
G  = spm_DEM_basis(x(1),v,P.p);
 
% cost function (C)
%--------------------------------------------------------------------------
C  = spm_DEM_basis(x(1),v,P.q);
 
% flow
%--------------------------------------------------------------------------
dt = 1/8;
f  = [x(2);
      G - x(2)*C]*dt;
 
 
% true scalar potential gradient (see spm_moutaincar_fx)
%--------------------------------------------------------------------------
% if x(1) < 0
%     G  = 2*x(1) + 1;
% else
%     xx = x(1)^2;
%     G  = (1 + 5*xx)^(-1/2) - 5*xx/(1 + 5*xx)^(3/2) + (x(1)/2)^4;
% end