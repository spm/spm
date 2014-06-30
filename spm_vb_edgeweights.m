function [edges,weights] = spm_vb_edgeweights(vxyz,img)
% Compute edge set and edge weights of a graph
% FORMAT [edges,weights]= spm_vb_edgeweights(vxyz,img)
% 
% vxyz     list of neighbouring voxels (see spm_vb_neighbors)
% img      image defined on the node set, e.g. wk_ols. The edge weights 
%          are uniform if this is not given, otherwise they are a function
%          of the distance in physical space and that between the image
%          at neighoring nodes

% edges    [Ne x 2] list of neighboring voxel indices
% weights  [Ne x 1] list of edge weights 
% Ne       number of edges (cardinality of edges set)
% N        number of nodes (cardinality of node set)
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Lee Harrison
% $Id: spm_vb_edgeweights.m 6079 2014-06-30 18:25:37Z spm $

N       = size(vxyz,1);
[r,c,v] = find(vxyz');
edges   = [c,v];
% undirected graph, so only need store upper [lower] triangle
i       = find(edges(:,2) > edges(:,1));
edges   = edges(i,:);
if nargin < 2,
    weights = ones(size(edges,1),1);
    return
else
    ka  = 16;
    M   = mean(img,2)*ones(1,N);
    C   = (1/N)*(img-M)*(img-M)';
    Hf  = inv(C);
    A   = spm_vb_incidence(edges,N);
    dB  = img*A'; % spatial gradients of ols estimates
    dg2 = sum((dB'*Hf).*dB',2); % squared norm of spatial gradient of regressors
    ds2 = 1 + dg2; % squared distance in space is 1 as use only nearest neighbors
    weights = exp(-ds2/ka); % edge weights
end
