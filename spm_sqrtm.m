function [K] = spm_sqrtm(V)
% Matrix square root (sqrtm) for sparse matrices
% FORMAT [K] = spm_sqrtm(V)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse
% matrices
 %___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sqrtm.m 2696 2009-02-05 20:29:48Z guillaume $

%--------------------------------------------------------------------------
n     = length(V);
[u s] = spm_svd(V,1e-16);
s     = sqrt(diag(s));
m     = length(s);
s     = sparse(1:m,1:m,s);
K     = u*s*u';
