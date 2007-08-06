function [D] = spm_speye(m,n,k)
% sparse leading diagonal matrix
% FORMAT [D] = spm_speye(m,n,k)
%
% returns an m x n matrix with ones along the k-th leading diagonal
%__________________________________________________________________________
 
% leading diagonal matrix
%--------------------------------------------------------------------------
D = spdiags(ones(m,1),k,m,n);
