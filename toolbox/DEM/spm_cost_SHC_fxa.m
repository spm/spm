function [f] = spm_cost_SHC_fxa(x,v,a,P)
% equations of motion for a foraging problem
% FORMAT [f] = spm_cost_SHC_fxa(x,v,a,P)
%
% x   - hidden states
% v   - exogenous inputs
% a   - action
% P   - parameters for mountain car
%
% returns f = dx/dt (see spm_cost_SHC_fx)
% These equations of motion model dissipative flow x.x and x.v on a flat 
% potential and increases in physiological states x.q as radial basis 
% functions of secrete locations. The agent has to discover these 
% locations % using an appropriate policy. This generative process would 
% also substitute for Morris water-maze simulations or unbounded saccades.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id:  $
 
% location and radius of attractors A (only A.q attractors deliver reward)
%--------------------------------------------------------------------------
global A
 
% physical flow
%--------------------------------------------------------------------------
f   = x;
f.x = x.v;
f.v = a - x.x/8 - x.v/4;
 
% physiological flow
%--------------------------------------------------------------------------
for i = 1:length(f.q)
    f.q(i) = exp(-sum((x.x - A.x(:,A.q(i))).^2)/(2*A.d^2)) - x.q(i)/2;
end
 
% flow
%--------------------------------------------------------------------------
dt  = 1/8;
f.x = f.x*dt;
f.v = f.v*dt;
f.q = f.q*dt;
