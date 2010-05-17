function [f]= spm_fx_dem_write(x,v,P)
% returns the flow for a two-joint arm (writing with SHC)
% FORMAT [f]= spm_fx_dem_write(x,v,P)
%
%   x.x(1) - joint angle
%   x.x(2) - joint angle
%   x.x(3) - angular velocity
%   x.x(4) - angular velocity
%
%   x.a(1) - attraction (location 1)
%   x.a(2) - attraction (location 2)
%   x.a(3) - attraction (location 3)
%    ...
%
% v    - hidden states
%   v(1) - not used
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_write.m 3893 2010-05-17 18:28:52Z karl $


% diameter of radial basis function (for autovitiation)
%--------------------------------------------------------------------------
d  = 1/8;

% joint location
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x.x);
X  = J{1} + J{2};

% motion of physical states
%==========================================================================

% desired location (P(:,i)) is determined by the attractor state x.a
%--------------------------------------------------------------------------
[m,i]  = max(x.a);
f.x    = spm_fx_dem_reach(x.x,P(:,i),P);
 
% motion of attractor states (using basis functions of position)
%==========================================================================
b      = zeros(size(x.a));
b(i,1) = norm(X - P(:,i)) < d;
f.a    = (-b*4 - sum(x.a))/8;
 



