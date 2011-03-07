function [f]= spm_fx_adem_cue(x,v,a,P)
% returns the flow for cued response (with action)
% FORMAT [f]= spm_fx_adem_cue(x,v,a,P)
%
% x    - hidden states:
%   x.o  - motor angle
%
% v    - hidden causes
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_adem_cue.m 4230 2011-03-07 20:58:38Z karl $
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;

% motion of oculomotor angles (driven by bounded action) with decay to 0
%==========================================================================
f.o  = tanh(a) - x.o/8;

% motion of target contrast
%==========================================================================
f.a  = v - x.a;