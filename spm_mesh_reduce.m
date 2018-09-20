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
% Wrapper around an implementation by Sven Forstmann, MIT licence:
% https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_reduce.m 7421 2018-09-20 10:58:01Z guillaume $


%-This is merely the help file for the compiled routine
% error('spm_mesh_reduce.cpp not compiled - see Makefile')

M = reducepatch(M,t);
