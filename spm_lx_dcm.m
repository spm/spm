function [y] = spm_lx_dcm(x,u,P)
% simulated BOLD response to input
% FORMAT [y] = spm_lx_dcm(x,u,P)
% y    - BOLD response (%)
%
% x    - state vector     (see spm_fx_HRF)
% P    - Parameter vector (see spm_fx_HRF)
%___________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%___________________________________________________________________________
% Karl Friston $Id$

% resting venous volume
%---------------------------------------------------------------------------
V0       = 100*0.02;
E0       = 0.34; 

% coeficients for BOLD signal
%---------------------------------------------------------------------------
k1       = 7*E0;
k2       = 2;
k3       = 2*E0 - 0.2;

% BOLD signal
%---------------------------------------------------------------------------
x        = reshape(x,length(x)/5,5);
x(:,3:5) = x(:,3:5) + 1;
y        = V0*(k1*(1 - x(:,5)) + k2*(1 - x(:,5)./x(:,4)) + k3*(1 - x(:,4)));
