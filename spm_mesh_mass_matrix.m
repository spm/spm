function A = spm_mesh_mass_matrix(M)
% Compute the mass matrix of a triangle mesh
% M        - patch structure: vertices and faces must be mx3 and nx3 arrays
%
% A        - Mass matrix
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_mass_matrix.m 8155 2021-09-26 16:29:44Z guillaume $


A = spm_mesh_area(M,'vertex');
A = spdiags(A,0,size(A,1),size(A,1));
