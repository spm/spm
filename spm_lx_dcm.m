function [y] = spm_lx_dcm(x,u,P)
% simulated BOLD response to input
% FORMAT [y] = spm_lx_dcm(x,u,P)
% y    - BOLD response (%)
% x    - state vector     (see spm_fx_HRF)
% P    - Parameter vector (see spm_fx_HRF)
%___________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_lx_dcm.m 294 2005-11-13 13:55:40Z klaas $


% resting venous volume
%---------------------------------------------------------------------------
V0       = 100*0.02;

% coeficients for BOLD signal
%---------------------------------------------------------------------------
% get hemodynamic parameters
m         = max(size(u));           % number of inputs
n         = max(size(x)/5);         % number of regions
[A B C H] = spm_dcm_reshape(P,m,n);
% estimated area-specific resting oxygen extraction fractions
E0        = H(:,5);
% compute coefficients k1,k2,k3 according to Buxton et al. (1998)
k1        = 7*E0;
k2        = 2 * ones(n,1);
k3        = 2*E0 - 0.2;

% BOLD signal
%---------------------------------------------------------------------------
x         = reshape(x,length(x)/5,5);
x(:,3:5)  = x(:,3:5) + 1;
y         = V0*(k1.*(1 - x(:,5)) + k2.*(1 - x(:,5)./x(:,4)) + k3.*(1 - x(:,4)));

return