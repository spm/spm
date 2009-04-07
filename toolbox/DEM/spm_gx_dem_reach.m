function [g]= spm_gx_dem_reach(x,v,P)
% returns the prediction for a two-joint arm
% FORMAT [g]= spm_gx_dem_reach(x,v,P)
%
% x    - hidden states
% v    - causal states
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_dem_reach.m 3054 2009-04-07 19:22:49Z karl $

% evaluate positions
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x);

% stretch (positional) and visual (positional) information
%==========================================================================
g  = [x(1:2); v; J{1} + J{2}];

