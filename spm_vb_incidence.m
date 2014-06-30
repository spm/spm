function A = spm_vb_incidence(edges,N)
% Edge-node incidence matrix of a graph
% FORMAT A = spm_vb_incidence(edges,N)
% 
% edges    [Ne x 2] list of neighboring voxel indices
% N        number of nodes (cardinality of node set)
%
% Ne       number of edges (cardinality of edge set)
% A        [Ne x N] matrix - is the discrete analogue of the grad operator
% A(ij,k)  +1 if i=k, -1 if j=k, 0 otherwise, where ij is the edge 
% connecting nodes i and j, and k is in node set
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Lee Harrison
% $Id: spm_vb_incidence.m 6079 2014-06-30 18:25:37Z spm $

% Number of edges
Ne = size(edges,1);

% Number of nodes (if N is not specified)
if nargin < 2
    N = max(edges(:));
end

% Edge-node incidence matrix
A   = sparse([1:Ne,1:Ne],[edges(:,1),edges(:,2)],...
        [ones(Ne,1),-ones(Ne,1)],Ne,N);
