
% returns voxel values from a memory mapped image - a compiled routine
% FORMAT X = spm_sample_vol(V,x,y,z,hold);
% V      -  is a memory mapped image volume
% x      -  matrix of x coordinates {pixels}
% y      -  matrix of y coordinates {pixels}
% z      -  matrix of z coordinates {pixels}
% hold   -  sets the interpolation method for the resampling.
%           0       Zero order hold (nearest neighbour).
%           1       First order hold (trilinear interpolation).
%           3-127   Sinc interpolation using different numbers of neighbours.
%                   Anything above about 8 is serious overkill
% X      -  output image
%____________________________________________________________________________
%
% spm_sample_vol will return the voxel values from a memory mapped volume
% indicated by V at coordinates x,y,z.  Values from coordinates outside the
% image are set to zero. x, y and z must be matrices of the same dimensions
%
% see also spm_slice_vol.m
%
%__________________________________________________________________________
% %W% %E%
