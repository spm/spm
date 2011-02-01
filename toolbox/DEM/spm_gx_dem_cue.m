function [g]= spm_gx_dem_cue(x,v,P)
% returns the prediction for cued responses (proprioception and vision)
% FORMAT [g]= spm_gx_dem_cue(x,v,P)
%
% x    - hidden states:
%   x.o  - oculomotor angle (proprioceptive)
%   x.x  - target locations (visual) - extrinsic coordinates (Cartesian)
%   x.a  - target contrast (attractiveness)
%
% v    - hidden causes
%
% P    - parameters
%
% g    - sensations:
%   g.o  - oculomotor angle (proprioception)
%   g.p  - target locations (visual) - intrinsic coordinates (polar)
%   g.c  - target contrast
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_dem_cue.m 4187 2011-02-01 20:13:57Z karl $
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = atan(x.x) - x.o*ones(1,size(x.x,2));
g.c = x.a;

