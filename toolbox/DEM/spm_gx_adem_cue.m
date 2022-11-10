function [g]= spm_gx_adem_cue(x,v,a,P)
% returns the prediction for cued responses (proprioception and vision)
% FORMAT [g]= spm_gx_adem_cue(x,v,a,P)
%
% x    - hidden states:
%   x.o  - intrinsic motor state (proprioceptive)
%
% v    - hidden causes
%
% P    - target locations (visual) - extrinsic coordinates (Cartesian)
%
% g    - sensations:
%   g.o  - motor angle (proprioception)
%   g.p  - finger location (visual)
%   g.c  - target contrast (visual)
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = tan(x.o);
g.c = x.a*4;

