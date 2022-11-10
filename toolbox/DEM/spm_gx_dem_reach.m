function [g]= spm_gx_dem_reach(x,v,P)
% returns the prediction for a two-joint arm
% FORMAT [g]= spm_gx_dem_reach(x,v,P)
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
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% evaluate positions
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x);

% stretch (angular) and visual (positional) information (target & arm)
%==========================================================================
g  = [x(1:2); v; J{1} + J{2}];

