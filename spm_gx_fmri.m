function [g,dgdx] = spm_gx_fmri(x,u,P,M)
% Simulated BOLD response to input
% FORMAT [g,dgdx] = spm_gx_fmri(x,u,P,M)
% g          - BOLD response (%)
% x          - state vector     (see spm_fx_fmri)
% P          - Parameter vector (see spm_fx_fmri)
% M          - model specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_fmri.m 6262 2014-11-17 13:47:56Z karl $
 
 
% Biophysical constants for 1.5T
%==========================================================================
 
% time to echo (TE) (default 0.04 sec)
%--------------------------------------------------------------------------
TE  = 0.04;
 
% resting venous volume (%)
%--------------------------------------------------------------------------
V0  = 4;

% estimated region-specific ratios of intra- to extra-vascular signal 
%--------------------------------------------------------------------------
ep  = exp(P.epsilon);
 
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
%--------------------------------------------------------------------------
r0  = 25;
 
% frequency offset at the outer surface of magnetized vessels (Hz)
%--------------------------------------------------------------------------
nu0 = 40.3; 
 
% resting oxygen extraction fraction
%--------------------------------------------------------------------------
E0  = 0.4;
 
%-Coefficients in BOLD signal model
%==========================================================================
k1  = 4.3*nu0*E0*TE;
k2  = ep*r0*E0*TE;
k3  = 1 - ep;
 
%-Output equation of BOLD signal model
%==========================================================================
v   = exp(x(:,4));
q   = exp(x(:,5));
g   = V0*(k1 - k1.*q + k2 - k2.*q./v + k3 - k3.*v);

if nargout == 1, return, end


%-derivative dgdx
%==========================================================================
[n m]      = size(x);
dgdx       = cell(1,m);
[dgdx{:}]  = deal(sparse(n,n));
dgdx{1,4}  = diag(-V0*(k3.*v - k2.*q./v));
dgdx{1,5}  = diag(-V0*(k1.*q + k2.*q./v));
dgdx       = spm_cat(dgdx);
