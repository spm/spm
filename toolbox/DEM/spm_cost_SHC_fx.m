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
% $Id: karl $
 
 
 
% location and radius of attractors A
%--------------------------------------------------------------------------
global A;
 
% gradient of Hamiltonian (G) is determined by the attractor state x.a
%--------------------------------------------------------------------------
f   = x;
n   = length(x.a);
G   = (x.x*ones(1,n) - A.x)*spm_softmax(x.a,2);
 
% motion of physical states
%--------------------------------------------------------------------------
f.x = x.v;
f.v = -G*8 - x.v*4;
 
% motion of physiological states (using basis functions of position)
%--------------------------------------------------------------------------
for i = 1:n
    b(i,1) = exp(-sum((x.x - A.x(:,i)).^2)/(2*A.d^2));
end
 
f.q = P'*b - x.q;
f.a = 1/2 + P*(x.q < A.u)*4 - b*8;
 
% flow
%--------------------------------------------------------------------------
dt  = 1/8;
f.x = f.x*dt;
f.v = f.v*dt;
f.q = f.q*dt;
f.a = f.a*dt;
