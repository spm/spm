function [f]= spm_fx_adem_salience(x,v,a,P)
% returns the flow for oculomotor search
% FORMAT [f]= spm_fx_adem_salience(x,v,a,P)
%
% x    - hidden states:
%   x(1) - oculomotor angle
%   x(2) - oculomotor angle
%
% v    - hidden cause
% P    - parameters
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 

% motion of oculomotor angles (driven by unbounded action)
%==========================================================================
f  = a - x/16;
