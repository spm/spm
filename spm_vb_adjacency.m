function W = spm_vb_adjacency(edges,weights,N)
% (Weighted) adjacency (or weight) matrix of a graph
% FORMAT W = spm_vb_adjacency(edges,weights,N)
%
% edges    [Nedges x 2] list of neighboring voxel indices
% weights  [Nedges x 1] list of edge weights (unity of not specified)
% N        number of nodes (cardinality of node set)
%
% W        [N x N] matrix of (weighted) edges
% Wij      edge weight between nodes i and j if they are neighbors, otherwise 0
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging
 
% Lee Harrison
% $Id: spm_vb_adjacency.m 6079 2014-06-30 18:25:37Z spm $

% Number of edges
Ne = size(edges,1);

% Uniform weights if not specified
if nargin < 2
    weights = ones(Ne,1);
end

% Number of nodes (if N is not specified)
if nargin < 3
    N = max(edges(:));
end

% (Weighted) adjacency matrix
W  = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],...
            [weights;weights],N,N);
