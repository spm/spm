function [I, t] = spm_mesh_ray_triangle(V, R)
% Compute the intersection of ray(s) and triangle(s)
% FORMAT [I, t] = spm_mesh_ray_triangle(V, R)
% V   - double precision matrix of triangles [v1;v2;v3] [9xNV]
% R   - double precision matrix of rays [origin;direction] [6xNR] 
%
% I   - face indices of intersection hit [NRxN]
% t   - distance to hit triangles [NRxN]
%__________________________________________________________________________
%
% This function implements the Moller-Trumbore ray-triangle intersection
% algorithm:
% "Fast, Minimum Storage Ray-Triangle Intersection". Tomas Moller and Ben
% Trumbore (1997). Journal of Graphics Tools. 2: 21-28.
% https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
% https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/
%__________________________________________________________________________
%
% This is a low-level function. Please use spm_mesh_ray_intersect.m
% instead.
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%-This is merely the help file for the compiled routine
error('spm_mesh_ray_triangle.c not compiled - see Makefile')
