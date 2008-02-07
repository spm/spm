function [K] = spm_sqrtm(V)
% matrix sqrt for sparse matrices
% FORMAT [K] = spm_sqrtm(V)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse
% matrices
 %___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sqrtm.m 1143 2008-02-07 19:33:33Z spm $

%--------------------------------------------------------------------------
n     = length(V);
[u s] = spm_svd(V,1e-16);
s     = sqrt(diag(s));
m     = length(s);
s     = sparse(1:m,1:m,s);
K     = u*s*u';
