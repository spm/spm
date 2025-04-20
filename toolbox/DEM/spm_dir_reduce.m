function [R] = spm_dir_reduce(a)
% Compression of a (Dirichlet) probability tensor
% FORMAT [R] = spm_dir_reduce(a)
% a  - [cell array] of likelhoods
% R  - reduction matrix: a*R = reduced Dirichlet likelihood mapping
%
% This routine returns a reduction matrix (R) that compresses a joint
% probability distribution parameterised with Dirichlet counts. This can be
% regarded as a simple and efficient kind of Bayesian model reduction in
% which the mutual information or expected free energy is preserved (i.e.,
% with minimal information loss). The reduction effectively combines
% columns of the implicit likelihood that are similar — in the sense that
% there information distance is less than one natural unit [more precisely,
% sqrt(2)]. This uses an efficient unique operator, based upon thresholded
% information distance.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
D       = spm_information_distance(a);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < sqrt(2),'rows','stable');

% restriction matrix for reduction of likelihood
%--------------------------------------------------------------------------
Ns  = numel(i);
R   = sparse(1:numel(j),j,1,numel(j),Ns);

% i.e.,
%--------------------------------------------------------------------------
% for i = 1:numel(j)
%     for j = 1:max(j)
%         if k(i) == j
%             R(i,j) = 1;
%         else
%             R(i,j) = 0;
%         end
%     end
% end



