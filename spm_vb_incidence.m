function A = spm_vb_incidence(edges,N)
% Edge-node incidence matrix of a graph
% FORMAT W = spm_vb_incidence(edges,N)
% 
% edges    [Ne x 2] list of neighboring voxel indices
% N        number of nodes (cardinality of node set)
%
% Ne       number of edges (cardinality of edge set)
% A        [Ne x N] matrix - is the discrete analogue of the grad operator
% A(ij,k)  +1 if i=k, -1 if j=k, 0 otherwise, where ij is the edge 
% connecting nodes i and j, and k is in node set
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Lee Harrison
% $Id: spm_vb_incidence.m 2921 2009-03-23 17:59:50Z guillaume $

% Number of edges
Ne = size(edges,1);

% Number of nodes (if N is not specified)
if nargin < 2
    N = max(edges(:));
end

% Edge-node incidence matrix
A   = sparse([1:Ne,1:Ne],[edges(:,1),edges(:,2)],...
        [ones(Ne,1),-ones(Ne,1)],Ne,N);