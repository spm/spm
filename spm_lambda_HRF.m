function [y] = spm_lambda_HRF(x,P)
% simulated BOLD response to input
% FORMAT [y] = spm_lambda_HRF(x,P)
% y    - BOLD response (%)
%
% x    - state vector     (see spm_fx_HRF)
% P    - Parameter vector (see spm_fx_HRF)
%___________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%___________________________________________________________________________
% %W% Karl Friston %E%

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
y     = 100*V0*(k1*(1 - x(4)) + k2*(1 - x(4)/x(3)) + k3*(1 - x(3)));
