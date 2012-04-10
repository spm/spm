function [K] = spm_sqrtm(V)
% Matrix square root for sparse symmetric positive semi-definite matrices
% FORMAT [K] = spm_sqrtm(V)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse
% symmetric positive semi-definite matrices.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sqrtm.m 4711 2012-04-10 13:20:39Z karl $

%--------------------------------------------------------------------------
[u s] = spm_svd(V,0,1e-16);
s     = sqrt(abs(diag(s)));
m     = length(s);
s     = sparse(1:m,1:m,s);
K     = u*s*u';
