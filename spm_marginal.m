function [Y] = spm_marginal(X)
% returns the marginal densities over a multidimensional array of probabilities
% FORMAT [Y] = spm_marginal(X)
% X  - numeric array of probabilities
%
% Y  - cell array marginals
%
% See also: spm_dot
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_marginal.m 7809 2020-03-31 11:55:09Z karl $

% use VECDIM in sum to evaluate marginals
%--------------------------------------------------------------------------
n     = ndims(X);
Y     = cell(n,1);
for i = 1:n
    j    = 1:n;
    j(i) = [];
    Y{i} = spm_vec(squeeze(sum(X,j)));
end
