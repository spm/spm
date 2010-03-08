function [f] = spm_cost_SHC_fx(x,v,P)
% equations of motion for foraging problem using SHCs
% problem
% FORMAT [f] = spm_cost_SHC_fx(x,v,P)
%
% x   - hidden states (x.x, x.v x.q and x.a)
% v   - exogenous inputs
% P   - parameters
%
% The parameters associate increases in some physiological states x.q with 
% positions in physical space, encoded by radial basis functions x.a
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cost_SHC_fx.m 3757 2010-03-08 11:41:53Z guillaume $
 
 
 
% location and radius of attractors A
%--------------------------------------------------------------------------
global A; X   = A.x;
 
% gradient of Hamiltonian (G) is determined by the attractor state x.a
%--------------------------------------------------------------------------
[m,i] = max(x.a);
G     = x.x - X(:,i);
 
% motion of physical states
%--------------------------------------------------------------------------
f   = x;
f.x = x.v;
f.v = -G*8 - x.v*4;
 
% motion of physiological states (using basis functions of position)
%--------------------------------------------------------------------------
for i = 1:size(X,2)
    b(i,1) = norm(x.x - X(:,i)) < A.d;
end
 
f.q = P'*b - x.q/2;
f.a = P*(x.q < A.u) - b*4 - sum(x.a);
 
% flow
%--------------------------------------------------------------------------
dt  = 1/8;
f.x = f.x*dt;
f.v = f.v*dt;
f.q = f.q*dt;
f.a = f.a*dt;
