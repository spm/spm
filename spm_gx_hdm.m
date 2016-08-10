function [y] = spm_gx_hdm(x,u,P,M)
% Simulated BOLD response to input.
% FORMAT [y] = spm_gx_hdm(x,u,P,M)
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
% $Id: spm_gx_hdm.m 6856 2016-08-10 17:55:05Z karl $


% biophysical constants for 1.5 T: 
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%--------------------------------------------------------------------------
if isstruct(P)
    H     = [0.64 0.32 2.00 0.32 0.4];
    for i = 1:numel(P.decay)
        H(6)   = P.epsilon;
        y(i,1) = spm_gx_hdm(x(i,:),u(i),H,M);
    end
    return
end


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
E0    = P(5); 

% region-specific ratios of intra- to extravascular components of
% the gradient echo signal (prior mean = 1, log-normally distributed 
% scaling factor)
%--------------------------------------------------------------------------
epsi  = exp(P(6));
 
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
v     = x(3);
q     = x(4);
y     = V0*(k1.*(1 - q) + k2.*(1 - (q./v)) + k3.*(1 - v)); 
