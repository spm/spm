function bb = spm_bb(P)
% Finds the 'bounding box' of an object in an image
% FORMAT bb = spm_bb(P)
% P     -       filename
% bb    -       [xmin ymin zmin
%                xmax ymax zmax] = dimensions of bounding box {voxels}
%___________________________________________________________________________
%
% spm_bb returns the bounding box containing an object in a volume
% image.   The box is defined by thresholding the integrated voxel
% values over each dimension at half the mean integrated counts
%
% see also spm_box.m
%
%__________________________________________________________________________
% %W% %E%

% get image dimensions and data type
%---------------------------------------------------------------------------
[DIM VOX SCALE TYPE] = spm_hread(P);

% get integrated values and set to a minimum of zero
%---------------------------------------------------------------------------
[x y z] = spm_box(P(P ~= ' '),DIM,TYPE);
x       = x - min(x);
y       = y - min(y);
u       = 2; 							% threshold

bb = [min(find(x >= mean(x)/u)) min(find(y >= mean(y)/u)) min(find(z >=mean(z)/u))
      max(find(x >= mean(x)/u)) max(find(y >= mean(y)/u)) max(find(z >=mean(z)/u))];

