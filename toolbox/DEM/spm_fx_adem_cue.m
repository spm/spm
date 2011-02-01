function [f]= spm_fx_adem_cue(x,v,a,P)
% returns the flow for cued response (with action)
% FORMAT [f]= spm_fx_adem_cue(x,v,a,P)
%
% x    - hidden states:
%   x.o  - oculomotor angle
%   x.a  - target contrast
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
% $Id: spm_fx_adem_cue.m 4187 2011-02-01 20:13:57Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;

% motion of oculomotor angles (driven by bounded action) with decay to 0
%==========================================================================
f.o  = tanh(a) - x.o/8;

% motion of target contrast
%==========================================================================
f.a  = spm_lotka_volterra(x.a,v);