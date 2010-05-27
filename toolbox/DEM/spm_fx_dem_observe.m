function [f]= spm_fx_dem_observe(x,v,P)
% returns the flow for a two-joint arm (writing with SHC)
% FORMAT [f]= spm_fx_dem_observe(x,v,P)
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
% $Id: spm_fx_dem_observe.m 3901 2010-05-27 16:14:36Z karl $

% motion of physical states
%==========================================================================

% desired location (v) is determined by the attractor state x.a
%--------------------------------------------------------------------------
n    = length(x.a);
p    = exp(2*x.a);
v    = P*p/sum(p);

% evaluate positions
%--------------------------------------------------------------------------
J  = spm_dem_reach_x2J(x.x);                  % joint location
F  = (v - J{1} - J{2})/2;                     % force


% flow
%==========================================================================
m    = [2  1]*2;                              % mass
k    = [2  1]*4;                              % viscosity
e    = 1/16;                                     % elasticity
O    = [0 -1 ;                                % orthogonal projector
        1  0];
f.x  = [x.x(3);
        x.x(4);
       (F'*J{2}*J{2}'*O*J{1} - k(1)*x.x(3) - (x.x(1) - pi/2)*e)/m(1);
       (           F'*O*J{2} - k(2)*x.x(4) - (x.x(2) - pi/2)*e)/m(2)];
 

% motion of (Lorenz) attractor states
%==========================================================================
T.f      = spm_speye(n,n,-1) - spm_speye(n,n,+1) + 2*spm_speye(n,n,0);
T.f(n,1) = -1;
T.f(1,n) = +1;
T.f      = T.f/2 - 1;
 
f.a      = spm_lotka_volterra(x.a,[],T)*.8;



