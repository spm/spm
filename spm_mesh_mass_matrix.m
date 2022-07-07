function A = spm_mesh_mass_matrix(M)
% Compute the mass matrix of a triangle mesh
% M        - patch structure: vertices and faces must be mx3 and nx3 arrays
%
% A        - Mass matrix
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


A = spm_mesh_area(M,'vertex');
A = spdiags(A,0,size(A,1),size(A,1));
