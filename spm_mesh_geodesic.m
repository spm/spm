function D = spm_mesh_geodesic(M,i,order)
% Approximate geodesic distances on a mesh using Dijkstra algorith
% FORMAT D = spm_mesh_geodesic(M,i,order)
% M        - a patch structure or an adjacency matrix
% i        - indices of starting points
% order    - type of distance for 1st order neighbours {0, [1]}
%
% D        - a [nx1] vector of geodesic distances from i
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_geodesic.m 4081 2010-10-07 14:04:44Z guillaume $

if nargin < 3, order = 1; end

[N, D] = spm_mesh_neighbours(M,order);

D      = spm_mesh_utils('dijkstra',N,D,i,Inf);
