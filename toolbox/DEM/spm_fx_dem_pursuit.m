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
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
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
t    = P*spm_softmax(x.a,1/2);
f.x  = (t - x.x)/2;

 
% motion of oculomotor angles (attracted to target)
%==========================================================================
t    = atan(x.x);
f.o  = (t - x.o);






