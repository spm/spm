function [g]= spm_gx_adem_write(x,v,a,P)
% returns the prediction for a two-joint arm (proprioception and vision)
% FORMAT [g]= spm_gx_adem_write(x,v,a,P)
%
% x    - hidden states:
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
% v    - causal states{
%   v(1) - exogenous force (x)
%   v(2) - exogenous force (y)
% a    - action
% P    - parameters
%
% g    - sensations:
%   g(1) - joint angle (proprioception)
%   g(2) - joint angle (proprioception)
%   g(3) - arm location (visual)
%   g(4) - arm location (visual)
% 
% As for spm_dem_reach but with no visual target
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_adem_write.m 3901 2010-05-27 16:14:36Z karl $

% evaluate positions
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x);

% stretch (angular) and visual (positional) information about motor plant
%==========================================================================
g  = [x; J{1}; J{1} + J{2}];

