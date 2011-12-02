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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_adem_salience.m 4580 2011-12-02 20:22:19Z karl $
 

% motion of oculomotor angles (driven by unbounded action)
%==========================================================================
f  = a - x/16;
