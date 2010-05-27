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
% P    - parameters (locations of point attratcors in state-space)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_write.m 3901 2010-05-27 16:14:36Z karl $


% diameter of radial basis function (for autovitiation)
%--------------------------------------------------------------------------
d  = 1/8;
n  = length(x.a);

% location in state-space
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x.x);
X  = J{1} + J{2};

% transition matrix
%--------------------------------------------------------------------------
T  = spm_speye(n,n,-1)/16 - spm_speye(n,n,0);


% motion of physical states
%==========================================================================

% desired location (P(:,i)) is determined by the attractor state x.a
%--------------------------------------------------------------------------
[m,i]  = max(x.a);
f.x    = spm_fx_dem_reach(x.x,[P(:,i); 1],P);

% motion of attractor states (using basis functions of position)
%==========================================================================
T      = T(:,i)*(norm(X - P(:,i)) < d);
f.a    = (T*6 - sum(x.a))/4;
 



