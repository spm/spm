function [A,i,j] = spm_dir_sort(A)
% sorts the rows and columns of a square matrix
% FORMAT [A,i,j] = spm_dir_sort(A)
%
% A    - matrix
% i,j  - indices
%
% Effectively, this reorders the rows and columns of A, so that the largest
% elements are along the leading diagonal of A(i,j)
%__________________________________________________________________________

% Karl Friston 
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% sort columns and then rows
%--------------------------------------------------------------------------
[m,i] = max(A);
[m,j] = sort(m,'descend');
[m,i] = max(A(:,j)');
[m,i] = sort(m,'descend');
A     = A(i,j);