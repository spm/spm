function [K] = spm_sqrtm(V)
% matrix sqrt for sparse matrices
% FORMAT [K] = spm_sqrtm(V)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse
% matrices
%__________________________________________________________________________
% %W% Karl Friston %E%
 
 
%--------------------------------------------------------------------------
n     = length(V);
[u s] = spm_svd(V,1e-16);
s     = sqrt(diag(s));
m     = length(s);
s     = sparse(1:m,1:m,s);
K     = u*s*u';
