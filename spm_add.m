
% add a series of images - a compiled routine
% FORMAT s = spm_add(V,Q,flags)
% V     - Vector of mapped volumes (from spm_map or spm_vol).
% Q     - Filename for averaged image
% flags - Flags can be:
%               'f' - writes floating point output image.
%               'm' - masks the mean to zero or NaN wherever
%                     a zero occurs in the input images.
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% spm_add computes a sum of a set of image volumes to produce an
% integral image that is written to a named file (Q).
%
% The image is written as signed short (16 bit) unless the `f' flag
% is specified.
%
% A mean can be effected by scaling the output image via it's
% scalefactor (see spm_mean for an example). A weighted sum can be
% effected by weighting the image scalefactors appropriately.
%
%_______________________________________________________________________
% %W% John Ashburner (from JB's original spm_mean.c) %E%
