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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_adem_cue.m 4230 2011-03-07 20:58:38Z karl $
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = tan(x.o);
g.c = x.a*4;

