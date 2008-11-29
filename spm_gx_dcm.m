function [y] = spm_gx_dcm(x,u,P,M)
% Simulated BOLD response to input.  This function implements the BOLD
% signal model described in Stephan et al. 2007, NeuroImage.
% FORMAT [y] = spm_gx_dcm(x,u,P,M)
% y    - BOLD response (%)
% x    - state vector     (see spm_fx_dcm)
% P    - Parameter vector (see spm_fx_dcm)
% M    - model  specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% Reference: 
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_dcm.m 2504 2008-11-29 15:53:11Z klaas $


% Biophysical constants for 1.5T
%--------------------------------------------------------------------------
% time to echo (TE)
try
    TE = M.TE;
catch
    TE = 0.04;
end
% resting venous volume
V0        = 100*0.04;
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
r0        = 25; % [Hz]
% frequency offset at the outer surface of magnetized vessels
nu0       = 40.3; % [Hz]


% Get estimates of hemodynamic parameters
%--------------------------------------------------------------------------
m         = max(size(u));           % number of inputs
n         = max(size(x)/5);         % number of regions
if ~M.nlDCM
    % bilinear DCM
    [A B C H]   = spm_dcm_reshape(P,m,n);
else
    % nonlinear DCM
    [A B C H D] = spm_dcm_reshape(P,m,n);
end

% estimated region-specific resting oxygen extraction fractions
E0        = H(:,5);

% estimated region-specific ratios of intra- to extravascular components of
% the gradient echo signal (prior mean = 1, log-normally distributed 
% scaling factor) 
epsilon   = exp(H(:,6));


% coefficients in BOLD signal model
%--------------------------------------------------------------------------
k1       = 4.3.*nu0.*E0.*M.TE;
k2       = epsilon.*r0.*E0.*M.TE;
k3       = 1 - epsilon;


% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x        = reshape(x,length(x)/5,5);
x(:,2:5) = exp(x(:,2:5)); 


% output equation of BOLD signal model
%--------------------------------------------------------------------------
v        = x(:,4);
q        = x(:,5);
% predicted BOLD response
y        = V0*(k1.*(1-q) + k2.*(1-(q./v)) + k3.*(1-v)); 

return
