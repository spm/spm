function [g] = spm_gx_adem_reach(x,v,a,P)
% returns the prediction for a two-joint arm (with action)
% FORMAT [g] = spm_gx_adem_reach(x,v,a,P)
%
% x    - hidden states
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
% v    - causal states
%   v(1) - target location (x)
%   v(2) - target location (y)
%   v(3) - force (cue strength)
% P    - parameters
% a    - action
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_adem_reach.m 3893 2010-05-17 18:28:52Z karl $

% stretch (angular) and visual (positional) information (target & arm)
%==========================================================================
g = spm_gx_dem_reach(x,v,P);