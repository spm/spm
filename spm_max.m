function varargout = spm_max(varargin)
% Sizes, maxima and locations of local excursion sets - a compiled routine
% FORMAT [N Z M A] = spm_max(X,L)
% X     - values of 3-D field
% L     - locations [x y x]' {in voxels}
% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in voxels}
% A     - region number
%_______________________________________________________________________
%
% spm_max characterizes a point list of voxel values (X) and their
% locations (L) in terms of edge, face and vertex connected subsets,
% returning a maxima- orientated list:  The value of the ith maximum is
% Z(i) and its location is given by M(:,i). A(i) identifies the ith
% maximum with a region. Region A(i) contains N(i) voxels.
%
% See also: spm_max.c and spm_clusters.m
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

%-This is merely the help file for the compiled routine
error('spm_max.c not compiled - see spm_MAKE.sh')
