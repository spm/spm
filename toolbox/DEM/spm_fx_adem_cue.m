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
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% intisaise flow (to ensure fields are aligned)
%--------------------------------------------------------------------------
f    = x;

% motion of oculomotor angles (driven by bounded action) with decay to 0
%==========================================================================
f.o  = tanh(a) - x.o/8;

% motion of target contrast
%==========================================================================
f.a  = v - x.a;