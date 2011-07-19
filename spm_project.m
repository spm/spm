function varargout = spm_project(varargin)
% forms maximium intensity projections - a compiled routine
% FORMAT spm_project(X,L,dims,[DXYZ,CXYZ])
% X -   a matrix of voxel values
% L -   a matrix of locations
% dims  -       assorted dimensions.
%               dims(1:3) - the sizes of the projected rectangles.
%               dims(4:5) - the dimensions of the mip image.
% optional
% DXYZ - length of the X,Y,Z axes of the mip sections (in mip pixels).
% CXYZ - offsets of the origin into the mip sections (in mip pixels).
%_______________________________________________________________________
%
% spm_project 'fills in' a matrix (SPM) in the workspace to create
% a maximum intensity projection according to a point list of voxel
% values (V) and their locations (L) in the standard space described
% in the atlas of Talairach & Tournoux (1988) or another space defined by
% a customised mip template.
%
% see also spm_mip.m and spm_mip_ui.m
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_project.m 4395 2011-07-19 08:44:42Z volkmar $


%-This is merely the help file for the compiled routine
error('spm_project.c not compiled - see Makefile')
