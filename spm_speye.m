function [D] = spm_speye(m,n,k)
% sparse leading diagonal matrix
% FORMAT [D] = spm_speye(m,n,k)
%
% returns an m x n matrix with ones along the k-th leading diagonal
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_speye.m 1131 2008-02-06 11:17:09Z spm $


% default k = 0
%--------------------------------------------------------------------------
if nargin < 3, k = 0; end
if nargin < 2, n = m; end
 
% leading diagonal matrix
%--------------------------------------------------------------------------
D = spdiags(ones(m,1),k,m,n);
