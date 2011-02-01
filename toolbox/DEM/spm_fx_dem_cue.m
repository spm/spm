function [f]= spm_fx_dem_cue(x,v,P)
% returns the flow for cued response
% FORMAT [f]= spm_fx_dem_cue(x,v,P)
%
% x    - hidden states:
%   x.o  - oculomotor angle
%   x.x  - target locations (visual) - extrinsic coordinates (Cartesian)
%   x.a  - target contrast (attractiveness)
%
% v    - hidden causes
%
% P    - parameters 
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_cue.m 4187 2011-02-01 20:13:57Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;

% motion of oculomotor angles (attracted to target)
%==========================================================================

% target location is determined by the attractor state softmax(x.a)
%--------------------------------------------------------------------------
t    = x.x*spm_softmax(x.a,2);
f.o  = (atan(t) - x.o)/2;

% motion of location states
%==========================================================================
f.x  = sparse(size(x.x,1),size(x.x,2));

% motion of attractor states
%==========================================================================
f.a  = spm_lotka_volterra(x.a,v(1));



