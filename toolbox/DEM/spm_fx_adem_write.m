function [f]= spm_fx_adem_write(x,v,a,P)
% returns the flow for a two-joint arm (with action)
% FORMAT [f]= spm_fx_adem_reach(x,v,a,P)
%
% x    - hidden states
%   x(1) - joint angle
%   x(2) - joint angle
%   x(3) - angular velocity
%   x(4) - angular velocity
% v    - exogenous forces (x,y)
% a    - action (forces) (x,y)
% P    - parameters
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% evaluate positions
%--------------------------------------------------------------------------
m    = [2  1]*2;                                 % mass
k    = [2  1]*4;                                 % viscosity

% flow
%==========================================================================
f    = [x(3);
        x(4);
      (-k(1)*x(3) + a(1) + v(1) - (x(1) - pi/2)/4)/m(1);
      (-k(2)*x(4) + a(2) + v(2) - (x(2) - pi/2)/4)/m(2)];
