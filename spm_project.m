
% forms maximium intensity projections - a compiled routine
% FORMAT spm_project(X,L,SPM,V)
% X	-	a matrix of voxel values
% L	- 	a matrix of locations in Talairach et Tournoux (1988) space
% SPM	-	matrix for maximum intensity projection
% V     -       {1 x 6} vector of image and voxel sizes [DIM VOX]
%____________________________________________________________________________
%
% spm_project 'fills in' a matrix (SPM) in the workspace to create
% a maximum intensity projection according to a point list of voxel
% values (X) and their locations (L) in the standard space described
% in the atlas of Talairach & Tournoux (1988).
%
% see also spm_mip.m
%
%__________________________________________________________________________
% %W% %E%
