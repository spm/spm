function [f] = spm_cost_fxa(x,v,a,P)
% equations of motion for a foraging problem
% problem
% FORMAT [f] = spm_cost_fxa(x,v,a,P)
%
% x   - hidden states
% v   - exogenous inputs
% a   - action
% P   - parameters for mountain car
%
% returns f = dx/dt 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cost_fxa.m 3757 2010-03-08 11:41:53Z guillaume $

% physical flow
%--------------------------------------------------------------------------
dt  = 1/8;
f   = x;
f.x = x.v;
f.v = a - x.v/8;

% physiological flow
%--------------------------------------------------------------------------
d     = 1/2;
for i = 1:size(P.p,2)
    n      = norm(x.x - P.p(:,i),'fro') < d;
    f.q(i) = n - x.q(i)/2;
    f.p(i) = x.q(i) - x.p(i)/8;
end

f.x = f.x*dt;
f.v = f.v*dt;
f.q = f.q*dt;
f.p = f.p*dt;


