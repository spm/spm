function varargout = spm_project(varargin)
% Form maximum intensity projections - a compiled routine
% FORMAT SPM = spm_project(X,L,dims,[DXYZ,CXYZ])
% X      -   a matrix of voxel values
% L      -   a matrix of locations
% dims   -   assorted dimensions.
%               dims(1:3) - the sizes of the projected rectangles.
%               dims(4:5) - the dimensions of the mip image.
% Optional:
% DXYZ   - length of the X,Y,Z axes of the mip sections (in mip pixels).
% CXYZ   - offsets of the origin into the mip sections (in mip pixels).
%__________________________________________________________________________
%
% spm_project 'fills in' a matrix (SPM) to create a maximum intensity
% projection according to a point list of voxel values (V) and their
% locations (L) in the standard space described in the atlas of Talairach &
% Tournoux (1988) or another space defined by a customised mip template.
%
% See also: spm_mip.m and spm_mip_ui.m
%__________________________________________________________________________
% Copyright (C) 1994-2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_project.m 5160 2012-12-21 16:58:38Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_project.c not compiled - see Makefile')
