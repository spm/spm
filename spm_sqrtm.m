function [K] = spm_sqrtm(V)
% Matrix square root for sparse symmetric positive semi-definite matrices
% FORMAT [K] = spm_sqrtm(V)
%
% This routine covers and extends sqrtm functionality by using a
% computationally expedient approximation that can handle sparse symmetric
% positive semi-definite matrices.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
[u,s] = spm_svd(V,0);
s     = sqrt(abs(diag(s)));
m     = length(s);
s     = sparse(1:m,1:m,s);
K     = u*s*u';
