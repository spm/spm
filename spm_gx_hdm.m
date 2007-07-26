function [y] = spm_gx_hdm(x,u,P,M)
% simulated BOLD response to input
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% y    - BOLD response (%) (y)
%
% x    - state vector      (see spm_fx_hdm)
% P    - Parameter vector  (see spm_fx_hdm)
%__________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_gx_hdm.m 868 2007-07-26 17:55:53Z karl $
 
 
% resting venous volume
%--------------------------------------------------------------------------
V0   = 0.04;
E0   = P(5); 
 
% coefficients for BOLD signal
%--------------------------------------------------------------------------
k1   = 7*E0;
k2   = 2;
k3   = 2*E0 - 0.2;
 
% exponentiation
%--------------------------------------------------------------------------
x    = exp(x);
 
% BOLD signal
%--------------------------------------------------------------------------
y(1) = 100*V0*(k1*(1 - x(4)) + k2*(1 - x(4)/x(3)) + k3*(1 - x(3)));
