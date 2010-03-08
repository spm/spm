function [f] = spm_cost_fx(x,v,P)
% equations of motion for foragaing problem
% problem
% FORMAT [f] = spm_cost_fx(x,v,P)
%
% x   - hidden states
% v   - exogenous inputs
% P.p - parameters for gradient function:     G(x(1),P.p)
%
% returns f = dx/dt = f  = [x(2);
%                           G - x(2)*C]*dt;
%
% where C determines divergence of flow x(2) at any position x(1).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cost_fx.m 3757 2010-03-08 11:41:53Z guillaume $
 

% dominant physiological state
%--------------------------------------------------------------------------
f     = x;
[p i] = min([x.p; 1]);

if i > length(x.p)
    
    A = [0;0];
    C = -1;
    
else
    
    A = P.p(:,i);
    q = x.q(i);
    
    % cost function
    %----------------------------------------------------------------------
    u = 1/16;
    C = (q < u)*1/8 - (q > u)*8;
    
end

% gradient of Hamiltonian (G)
%--------------------------------------------------------------------------
k   = 2;
G   = (x.x - A)*k;

% flow
%--------------------------------------------------------------------------
f.x = x.v;
f.v = x.v*C - G;

% physiological flow
%--------------------------------------------------------------------------
d   = 1/2;
for i = 1:size(P.p,2)
    f.q(i) = (norm(x.x - P.p(:,i),'fro') < d) - x.q(i);
    f.p(i) = x.q(i) - x.p(i)/8;
end

dt  = 1/8;
f.x = f.x*dt;
f.v = f.v*dt;
f.q = f.q*dt;
f.p = f.p*dt;


