function [i,j] = spm_unique(O)
% Information geometry of a likelihood mapping or [probabilitic] outcome
% FORMAT [i,j] = spm_unique(O)
% O{n,m}    - cell array of Dirichlet counts or outcomes
%
% returns indices (i,j) where U = O(i) and O = U(j)
% where U = unique(O{:,u})
%
% This routine uses the information geometry inherent in a likelihood
% mapping. It computes a distance matrix based upon the information length
% between columns of likelihood tensors and finds unique columns.
%__________________________________________________________________________
% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% Fast approximation by simply identifying unique locations on a
% multinomial statistical manifold, using Matlab's unique operator.
%--------------------------------------------------------------------------
O       = spm_dir_norm(O);
[~,i,j] = unique(logical(spm_cat(O)'),'rows','stable');


return

% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
D       = spm_information_distance(O);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < 2,'rows','stable');

