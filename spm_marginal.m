function [Y] = spm_marginal(X)
% Marginal densities over a multidimensional array of probabilities
% FORMAT [Y] = spm_marginal(X)
% X  - numeric array of probabilities
%
% Y  - cell array of marginals
%
% See also: spm_dot
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


% evaluate marginals
%--------------------------------------------------------------------------
n     = ndims(X);
Y     = cell(n,1);
for i = 1:n
    j    = 1:n;
    j(i) = [];
    Y{i} = reshape(spm_sum(X,j),[],1);
end
