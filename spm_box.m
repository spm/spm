
% integrates a volume image over x, y and z - a compiled routine
% FORMAT [X Y Z] = spm_box(P,DIM,TYPE)
% P	-	filename (unsigned char)
% DIM	-	[x y z] - image size {voxels}
% TYPE  -       data type (see spm_type.m & volume.h)
% X,Y,Z	-	integrated 1-dimensional images
%____________________________________________________________________________
%
% spm_box simply integrates (sums) all voxel values over x, y and z to
% give three vectors of integrated voxel values.  It is used primarily
% to determine the bounding box for the 'object' in image space.
%
% see also spm_bb.m
%
%__________________________________________________________________________
% %W% %E%
