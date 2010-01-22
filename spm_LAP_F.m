function [L] = spm_LAP_F(q,qu,qp,qh,pu,pp,ph,M)
% returns the Gibbs energy (L) as a function of contitional means
% FORMAT [L] = spm_LAP_F(q,qu,qp,qh,M)
%
%     q.x: {nx1 cell}
%     q.v: {dx1 cell}
%     q.p: {mx1 cell}
%     q.h: {mx1 cell}
%     q.g: {mx1 cell}
%
% for an m-level hierarchy
% See spm_LAP
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_LAP_F.m 3694 2010-01-22 14:16:51Z karl $

% Place conditional expecations in qu, qp, qh
%--------------------------------------------------------------------------
n  = length(q.x);
d  = length(q.v);

qu.x(1:n) = q.x;
qu.v(1:d) = q.v;
qp.p      = q.p;
qh.h      = q.h;
qh.g      = q.g;

% prediction errors
%--------------------------------------------------------------------------
E   = spm_DEM_eval(M,qu,qp);

Eu  = spm_vec(qu.x(1:n),qu.v(1:d));
Ep  = spm_vec(qp.p);
Eh  = spm_vec(qh.h,qh.g) - ph.h;

% precision and covariance matrices
%--------------------------------------------------------------------------
p   = spm_LAP_eval(M,qu,qh);
R   = spm_DEM_R(n,M(1).E.s);
iS  = spm_LAP_iS(p,R);

% Free-energy
%--------------------------------------------------------------------------
L   = E'*iS*E/2      - spm_logdet(iS)/2    + ...
      Eu'*pu.ic*Eu/2 - spm_logdet(pu.ic)/2 + ...
      Ep'*pp.ic*Ep/2 - spm_logdet(pp.ic)/2 + ...
      Eh'*ph.ic*Eh/2 - spm_logdet(ph.ic)/2;
