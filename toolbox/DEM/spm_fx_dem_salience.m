function [f]= spm_fx_dem_salience(x,v,P)
% returns the flow for visual search
% FORMAT [f]= spm_fx_dem_salience(x,v,P)
%
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes - attracting location
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_salience.m 4862 2012-08-24 19:21:27Z karl $
 
% intialise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;
 
% motion of oculomotor angles (attracted to target)
%==========================================================================
f.o  = (v - x.o)/4;

% motion of hypothesis states (with competition)
%==========================================================================
f.x = 1 - sum(exp(x.x)) - x.x/128;

