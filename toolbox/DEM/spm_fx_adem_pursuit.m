function [f]= spm_fx_adem_pursuit(x,v,a,P)
% returns the flow for occulomotor pursuit (with action)
% FORMAT [f]= spm_fx_adem_pursuit(x,v,a,P)
%
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.a(:) - attractor (SHC) states
%
% v    - hidden cause (speed)
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_adem_pursuit.m 4625 2012-01-24 20:53:10Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;
 
% motion of attractor states
%==========================================================================
f.a  = spm_lotka_volterra(x.a,v);
 
 
% motion of target states
%==========================================================================
 
% target location is determined by the attractor state x.a
%--------------------------------------------------------------------------
t    = P*spm_softmax(x.a,1/2);
f.x  = (t - x.x)/2;
 
 
% motion of oculomotor angles (driven by bounded action)
%==========================================================================
f.o  = tanh(a) - x.o/16;
