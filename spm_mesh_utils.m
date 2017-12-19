function varargout = spm_mesh_utils(action,varargin)
% A gateway function for surface mesh-related compiled algorithms 
%
% FORMAT [N, D] = spm_mesh_utils('neighbours',A)
% Return an array of first-order neighbours given an adjacency matrix
%
% FORMAT Fi = spm_mesh_utils('neighbouringfaces',F,i)
% Return the indices of the neighbouring triangles of a given triangle
% 
% FORMAT D = spm_mesh_utils('dijkstra',N,D,i,dmax)
% Compute geodesic distance on a triangular mesh using Dijkstra algorith
%
% FORMAT V = spm_mesh_utils('volume',M)
% Compute the volume of a closed surface mesh
%__________________________________________________________________________
% Copyright (C) 2010-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_utils.m 7240 2017-12-19 12:06:59Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_mesh_utils.c not compiled - see Makefile')
