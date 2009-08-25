function [f] = spm_mc_fx_5(x,v,P)
% equations of motion for the mountain car problem using basis functions
% problem
% FORMAT [f] = spm_mc_fx_5(x,v,P)
%
% x   - hidden states
% v   - exogenous inputs
% P.p - parameters for gradient function:     G(x(1),P.p)
% P.q - parameters for cost or loss-function: C(x(1),P.q)
%
% returns f = dx/dt = f  = [x(2);
%                           G - x(2)*C]*dt;
%
% where C determines divergence of flow x(2) at any position x(1).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mc_fx_5.m 3333 2009-08-25 16:12:44Z karl $
 
 
% gradient of Hamiltonian (G)
%--------------------------------------------------------------------------
G   = spm_DEM_basis(x.x,v,P);

% cost function (C)
%--------------------------------------------------------------------------
D   = tanh((1 - x.p)*8);
C   = 2 - 32*x.q;
C   = C*D - 1;

% flow
%--------------------------------------------------------------------------
dt  = 1/8;
f.x = [x.v; G - x.v*x.c]*dt;
f.c = [-C - x.c]*dt;
f.p = -x.p/64;
f.q = -x.q/64;
 
 

return

% NOTES for graphics
%--------------------------------------------------------------------------
x     = -2:1/64:2;
d     =  0:1/64:2;
for i = 1:length(x)
    for j = 1:length(d)
        D      = spm_phi((1 - d(j))*8);
        A      = 2 - 32*exp(-(x(i) - 1).^2*32);
        C(i,j) = A*D - 1;
    end
end

surf(d,x,C)
shading interp
xlabel('drive')
ylabel('position')

% true scalar potential gradient (see spm_moutaincar_fx)
%--------------------------------------------------------------------------
% if x(1) < 0
%     G  = 2*x(1) + 1;
% else
%     xx = x(1)^2;
%     G  = (1 + 5*xx)^(-1/2) - 5*xx/(1 + 5*xx)^(3/2) + (x(1)/2)^4;
% end