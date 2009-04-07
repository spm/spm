function [g] = spm_gx_adem_reach(x,v,a,P)
% returns the prediction for a two-joint arm (with action)
% FORMAT [g] = spm_gx_adem_reach(x,v,a,P)
%
% x    - hidden states
% v    - causal states
% a    - action
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_adem_reach.m 3054 2009-04-07 19:22:49Z karl $

% stretch (positional) and visual (positional) information
%==========================================================================
g = spm_gx_dem_reach(x,v,P);