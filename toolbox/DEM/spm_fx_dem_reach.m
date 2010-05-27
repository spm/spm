function [f]= spm_fx_dem_reach(x,v,P)
% returns the flow for a two-joint arm
% FORMAT [f]= spm_fx_dem_reach(x,v,P)
%
% x    - hidden states
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
% v    - causal states
%   v(1) - target location (x)
%   v(2) - target location (y)
%   v(3) - cue strength
% P    - parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_dem_reach.m 3901 2010-05-27 16:14:36Z karl $


% evaluate positions
%--------------------------------------------------------------------------
m  = [2  1]*2;                                % mass
k  = [2  1]*4;                                % friction
O  = [0 -1 ;                                  % orthogonal projector
      1  0];

T  = [v(1); v(2)];                            % target location
J  = spm_dem_reach_x2J(x);                    % joint location
F  = v(3)*(T - J{1} - J{2})*2;                  % force


% flow
%==========================================================================
f  = [x(3);
      x(4);
     (F'*J{2}*J{2}'*O*J{1} - k(1)*x(3))/m(1);
     (           F'*O*J{2} - k(2)*x(4))/m(2)];
