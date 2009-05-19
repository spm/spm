function C = spm_mesh_curvature(M)
% Compute a crude approximation of the curvature of a surface mesh
% FORMAT C = spm_mesh_curvature(M)
% M        - a patch structure
%
% C        - curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_curvature.m 3135 2009-05-19 14:49:42Z guillaume $

A = spm_mesh_adjacency(M);
A = sparse(1:size(M.vertices,1),1:size(M.vertices,1),1./sum(A,2)) * A;

C = (A-speye(size(A))) * double(M.vertices);
N = spm_mesh_normals(M);
C = sign(sum(N.*C,2)) .* sqrt(sum(C.*C,2));
