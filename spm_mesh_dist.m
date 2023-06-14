function D = spm_mesh_dist(M, XYZ)
% Compute signed distance to a triangle mesh
% FORMAT D = spm_mesh_dist(M, XYZ)
% M        - a patch structure with fields 'faces' and 'vertices'
% XYZ      - a n x 3 array of coordinates {mm}
%
% D        - a n x 1 vector of signed distances to the triangle mesh
%__________________________________________________________________________
%
% Based on C++ library:
% https://github.com/InteractiveComputerGraphics/TriangleMeshDistance
% Copyright (c) 2021 Jose Antonio Fernandez Fernandez, MIT license
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
error('spm_mesh_dist.cpp not compiled - see Makefile')
