function varargout = spm_project(varargin)
% forms maximium intensity projections - a compiled routine
% FORMAT spm_project(X,L,dims)
% X -   a matrix of voxel values
% L -   a matrix of locations in Talairach et Tournoux (1988) space
% dims  -       assorted dimensions.
%               dims(1:3) - the sizes of the projected rectangles.
%               dims(4:5) - the dimensions of the mip image.
%_______________________________________________________________________
%
% spm_project 'fills in' a matrix (SPM) in the workspace to create
% a maximum intensity projection according to a point list of voxel
% values (V) and their locations (L) in the standard space described
% in the atlas of Talairach & Tournoux (1988).
%
% see also spm_mip.m
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_project.m 1143 2008-02-07 19:33:33Z spm $


%-This is merely the help file for the compiled routine
error('spm_project.c not compiled - see Makefile')
