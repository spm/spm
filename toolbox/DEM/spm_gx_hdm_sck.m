function [y] = spm_gx_hdm_sck(x,u,P,M)
% Simulated BOLD response to input.
% FORMAT [y] = spm_gx_hdm_sck(x,u,P,M)
% y    - BOLD response (%)
% x    - state vector     (see spm_fx_fmri)
% P    - Parameter vector (see spm_fx_fmri)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_hdm_sck.m 6263 2014-11-17 13:48:36Z karl $


% biophysical constants for 1.5 T: 
%==========================================================================

% echo time (seconds)
%--------------------------------------------------------------------------
try
    TE = M(1).TE;
catch
    TE = 0.04;
end

% resting venous volume
%--------------------------------------------------------------------------
V0    = 100*0.08;                                

% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
%--------------------------------------------------------------------------
r0    = 25;

% frequency offset at the outer surface of magnetized vessels
%--------------------------------------------------------------------------
nu0   = 40.3;

% region-specific resting oxygen extraction fractions
%-------------------------------------------------------------------------- 
E0    = P(5,:); 

% region-specific ratios of intra- to extravascular components of
% the gradient echo signal (prior mean = 1, log-normally distributed 
% scaling factor)
%--------------------------------------------------------------------------
epsi  = exp(P(6,:));
 
% coefficients in BOLD signal model
%--------------------------------------------------------------------------
k1    = 4.3.*nu0.*E0.*TE;
k2    = epsi.*r0.*E0.*TE;
k3    = 1 - epsi;
 
% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x     = exp(x); 

% BOLD signal
%--------------------------------------------------------------------------
v     = x(3,:);
q     = x(4,:);
y     = V0.*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v)); 
