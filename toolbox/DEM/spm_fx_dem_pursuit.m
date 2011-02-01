function [f]= spm_fx_dem_pursuit(x,v,P)
% returns the flow for visual pursuit demo
% FORMAT [f]= spm_fx_dem_pursuit(x,v,P)
%
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.a(:) - attractor (SHC) states
%
% v    - hidden causes
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_pursuit.m 4187 2011-02-01 20:13:57Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;
 
% motion of attractor states
%==========================================================================
f.a  = spm_lotka_volterra(x.a,v);
 
 
% motion of target states
%==========================================================================
 
% target location is determined by the attractor state softmax(x.a)
%--------------------------------------------------------------------------
n    = length(x.a);
p    = exp(2*x.a);
t    = P*p/sum(p);
f.x  = (t - x.x)/8;

 
% motion of oculomotor angles (attracted to target)
%==========================================================================
t    = atan(x.x);
f.o  = (t - x.o)/2;






