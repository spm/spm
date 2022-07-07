function C = spm_mesh_curvature(M)
% Compute a crude approximation of the curvature of a surface mesh
% FORMAT C = spm_mesh_curvature(M)
% M        - a patch structure
%
% C        - curvature vector
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging


A = spm_mesh_adjacency(M);
A = sparse(1:size(M.vertices,1),1:size(M.vertices,1),1./sum(A,2)) * A;

C = (A-speye(size(A))) * double(M.vertices);
N = spm_mesh_normals(M);
C = sign(sum(N.*C,2)) .* sqrt(sum(C.*C,2));
