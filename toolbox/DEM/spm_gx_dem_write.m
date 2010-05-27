function [g]= spm_gx_dem_write(x,v,P)
% returns the prediction for a two-joint arm (writing example)
% FORMAT [g]= spm_gx_dem_write(x,v,P)
%
% x    - hidden states:
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
%
% v    - hidden causes
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
% $Id: spm_gx_dem_write.m 3901 2010-05-27 16:14:36Z karl $

% evaluate positions
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x.x);

% stretch (angular) and visual (positional) information about motor plant
%==========================================================================
g  = [x.x; J{1}; J{1} + J{2}];

