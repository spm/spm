
% integrates a series of images - a compiled routine
% FORMAT s = spm_mean(V,Q,flags)
% V     - Vector of mapped volumes (from spm_map or spm_vol).
% Q     - Filename for averaged image
% flags - Flags can be:
%               'f' - writes floating point output image.
%               'm' - masks the mean to zero or NaN wherever
%                     a zero occurs in the input images.
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% spm_mean computes a sum of a set of image volumes to produce an
% integral image that is written to a named file (Q).
% The image is written as signed short (16 bit) unless the `f' flag
% is specified.
%
% See also: spm_average - for details on how to usefully use this function
%
%_______________________________________________________________________
% %W% John Ashburner (highly modified from JB's original version) %E%
