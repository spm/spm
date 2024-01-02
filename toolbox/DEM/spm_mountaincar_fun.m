function [f] = spm_mountaincar_fun(P,G)
% [Cross-entropy] objective function for mountain car problem
% FORMAT [f] = spm_mountaincar_fun(P,G)
%
% P = spm_vec(P)
% P.a - 0th order coefficients of force
% P.b - 1st order coefficients of force
% P.c - 2nd order coefficients of force
% P.d - action efficacy
%
% G   - world model; including
%    G.fq : function fq(x) returning desired equilibrium density at x
%    G.X  : matrix of locations in x
%
% f   - KL divergence between actual and desired equilibrium densities
% x   - cell of grid point support 
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraße 38, 72076 T¨ubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
 
% place parameters in model
%--------------------------------------------------------------------------
G(1).pE = P;
 
% get equilibrium density
%--------------------------------------------------------------------------
[M0,q0] = spm_fp(G);
 
% desired equilibrium density
%--------------------------------------------------------------------------
q   = feval(G(1).fq,G(1).X);
  
% KL divergence or cross entropy
%--------------------------------------------------------------------------
q   = q/sum(q(:));
D   = q(:)./q0(:);
i   = find(D > exp(-16));
f   = q(i)'*log(D(i));
