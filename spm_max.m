function [N Z M A] = spm_max(X,L,VOX)
% Return sizes, maxima and locations of local excursion {X > u} sets 
% FORMAT [N Z M A] = spm_max(X,L,VOX)
% X     - values of 3-D field
% L     - locations [x y x]' {in mm}
% VOX   - voxel size {in mm}
% N     - size of region {in voxels)
% Z     - Z values of maxima
% M     - location of maxima {in mm}
% A     - region number
%____________________________________________________________________________
%
% spm_max characterizes a point list of voxel values (X) and their locations
% (L) in terms of edge, face and vertex connected subsets, returning a maxima-
% orientated list:  The value of the ith maximum is Z(i) and its location
% is given by M(:,i). A(i) identifies the ith maximum with a region. Region
% A(i) contains N(i) voxels.
%
% see also spm_max.c and spm_clusters.m
%
%__________________________________________________________________________
% %W% %E%
