function [g]= spm_gx_adem_cue(x,v,a,P)
% returns the prediction for cued responses (proprioception and vision)
% FORMAT [g]= spm_gx_adem_cue(x,v,a,P)
%
% x    - hidden states:
%   x.o  - oculomotor angle
%   x.a  - target contrast (attractiveness)
%
% v    - hidden causes
%
% P    - parameters 
%  P.x   - target locations - extrinsic coordinates (Cartesian)
%
% g    - sensations:
%   g.o  - oculomotor angle (proprioception)
%   g.p  - target locations (visual) - intrinsic coordinates (polar)
%   g.c  - target contrast
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_adem_cue.m 4187 2011-02-01 20:13:57Z karl $
 
% evaluate positions in intrinsic (polar) coordinates
%--------------------------------------------------------------------------
g.o = x.o;
g.p = atan(P.x) - x.o*ones(1,size(P.x,2));
g.c = x.a;

