function D = spm_mesh_dist(M, XYZ, S)
% Compute signed or unsigned distance from a point to a triangle mesh
% FORMAT D = spm_mesh_dist(M, XYZ, S)
% M        - a patch structure with fields 'faces' and 'vertices'
% XYZ      - a n x 3 array of points coordinates {mm}
% S        - logical scalar for signed or unsigned distances
%            [default: true, i.e. signed]
%
% D        - a n x 1 vector of signed or unsigned distances from XYZ to M
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
