function M = spm_mesh_reduce(M,t)
% Reduce the number of triangles in a mesh
% FORMAT M = spm_mesh_reduce(M,f)
% M        - a patch structure
% t        - desired number of triangles
%
% M        - reduced patch structure
%__________________________________________________________________________
%
% References:
%
% M. Garland and P. Heckbert. Surface Simplification Using Quadric Error
% Metrics. In Proceedings of SIGGRAPH 97.
% http://mgarland.org/files/papers/quadrics.pdf
% 
% M. Garland and P. Heckbert. Simplifying Surfaces with Color and Texture
% using Quadric Error Metrics. In Proceedings of IEEE Visualization 98. 
% http://mgarland.org/files/papers/quadric2.pdf
% 
% Wrapper around a C++ implementation by Sven Forstmann, MIT licence:
% https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
% ported to pure C by Chris Rorden, BSD 2-Clause License:
% https://github.com/neurolabusc/nii2mesh
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2018-2022 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
% error('spm_mesh_reduce.c not compiled - see Makefile')

M = reducepatch(M,t);
