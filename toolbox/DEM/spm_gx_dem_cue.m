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
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = tan(x.o);
g.c = exp(x.a/2);
