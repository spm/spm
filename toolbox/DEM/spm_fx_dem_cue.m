function [f]= spm_fx_dem_cue(x,v,P)
% returns the flow for cued response
% FORMAT [f]= spm_fx_dem_cue(x,v,P)
%
% x    - hidden states:
%   x.o  - intrinsic motor state (proprioceptive)
%   x.a  - target salience (attractiveness)
%
% v    - hidden causes
%
% P.x  - target locations (visual) - extrinsic coordinates (Cartesian)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_cue.m 4230 2011-03-07 20:58:38Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;

% motion of oculomotor angles (attracted to target)
%==========================================================================

% target location is determined by the attractor state softmax(x.a)
%--------------------------------------------------------------------------
t    = P.x*spm_softmax(x.a,P.s);
f.o  = (t - tan(x.o))/2;

% motion of attractor states
%==========================================================================
f.a  = spm_lotka_volterra(x.a,1/2*v(1));


