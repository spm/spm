function [g]= spm_gx_dem_cue(x,v,P)
% returns the prediction for cued responses (proprioception and vision)
% FORMAT [g]= spm_gx_dem_cue(x,v,P)
%
% x    - hidden states:
%   x.o  - intrinsic motor state (proprioceptive)
%   x.a  - target salience (attractiveness)
%
% v    - hidden causes
%
% P.x  - target locations (visual) - extrinsic coordinates (Cartesian)
%
% g    - sensations:
%   g.o  - motor angle (proprioception)
%   g.p  - finger locations (visual)
%   g.c  - target contrast  (visual)
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_dem_cue.m 4230 2011-03-07 20:58:38Z karl $
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = tan(x.o);
g.c = exp(x.a/2);
