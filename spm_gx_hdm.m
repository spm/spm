function [y] = spm_gx_hdm(x,u,P,M)
% simulated BOLD response to input
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
% y    - BOLD response (%) (y)
%
% x    - state vector      (see spm_fx_hdm)
% P    - Parameter vector  (see spm_fx_hdm)
%___________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_gx_hdm.m 740 2007-02-16 20:56:22Z karl $


% resting venous volume
%---------------------------------------------------------------------------
V0    = 0.02;
E0    = P(5); 

% coeficients for BOLD signal
%---------------------------------------------------------------------------
k1    = 7*E0;
k2    = 2;
k3    = 2*E0 - 0.2;

% BOLD signal
%---------------------------------------------------------------------------
y(1)  = 100*V0*(k1*(1 - x(4)) + k2*(1 - x(4)/x(3)) + k3*(1 - x(3)));

