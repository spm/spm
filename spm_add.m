
% add a series of images - a compiled routine
% FORMAT s = spm_add(VI,VO)
% VI    - Vector of mapped volumes (from spm_map or spm_vol).
% VO    - Description of output volume that gets passed to
%         spm_write_plane.m
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% spm_add computes a sum of a set of image volumes to produce an
% integral image that is written to a named file (VI.fname).
%
% A mean can be effected by modifying the scalefactors (and offsets) of
% VI (see spm_mean for an example). A weighted sum can be effected by
% using different weightings for image scalefactors.
%
% See also, spm_res2.m
%_______________________________________________________________________
% %W% John Ashburner (from JB's original spm_mean.c) %E%
