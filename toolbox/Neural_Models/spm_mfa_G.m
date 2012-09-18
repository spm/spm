function [G] = spm_mfa_G(M,x)
% creates a structure for a Gibb's ensemble
% FORMAT [G] = spm_mfa_G(M,x)
%--------------------------------------------------------------------------
% M   - model specification structure
% Required fields:
%    M.f   - dx/dt    = f(x,u,P)                {function string or m-file}
%    M.g   - y(t)     = l(x,P)                  {function string or m-file}
%    M.m   - m inputs
%    M.n   - n states
%    M.l   - l ouputs
%    M.x   - (n x 1) = x(0) = expansion point
%    M.W   - (n x n) - covariance matrix of deterministic noise
% x    -  {1 x d cell}    - range of state [d]-space
%
% G   - ensemble specification structure
% fields
%   G.M:  [1 x 1 struct]  - dynamic model structure
%   G.J0: [n x n double]  - Jacobian
%   G.J1: {1 x M.m cell}  - dJ0/du
%   G.L : [l x n double]  - d<y>/dp
%   G.u:  [n x m double]  - probability modes
%   G.v:  [m x n double]  - v*u = 1
%   G.X:  [n x d double]  - evaluation points of state [d]-space
%   G.x:  {1 x d cell}    - range of state [d]-space
%   G.p0: [n x 1 sparse]  - expansion point
%   G.q0: [n x 1 sparse]  - equilibrium density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mfa_G.m 4936 2012-09-18 19:47:55Z karl $
 
% Bilinear Mean field form
%--------------------------------------------------------------------------
[J0,J1,L,X,q0] = spm_mfa(M,x);
 
% find modes (remove unstable modes)
%--------------------------------------------------------------------------
if length(J0) < 1024
    [u s] = eig(full(J0));
else
    [u s] = eigs(J0,65,'LR');
end
r     = real(diag(s));
[i j] = sort(-r);
p0    = u(:,j(1));
u     = u(:,j((1:64) + 1));
v     = pinv(u);
 
% create structure for Gibb's ensemble
%--------------------------------------------------------------------------
G.M   = M;
G.J0  = J0;
G.J1  = J1;
G.L   = L;
G.u   = u;
G.v   = v;
G.x   = x;
G.X   = X;
G.p0  = p0/sum(p0);
G.q0  = q0/sum(q0);


