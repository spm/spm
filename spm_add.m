function varargout = spm_add(varargin)
% add a series of images - a compiled routine
% FORMAT s = spm_add(VI,VO)
% VI    - Vector of mapped volumes (from spm_map or spm_vol).
% VO    - Description of output volume that gets passed to
%         spm_write_plane.m
% flags - Flags can be:
%               'm' - masks the mean to zero or NaN wherever
%                     a zero occurs in the input images.
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% spm_add computes a sum of a set of image volumes to produce an
% integral image that is written to a named file (VI.fname).
%
% A mean can be effected by modifying the scalefactors (and offsets) of
% VI (see spm_mean_ui for an example). A weighted sum can be effected by
% using different weightings for image scalefactors.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_add.m 1143 2008-02-07 19:33:33Z spm $


%-This is merely the help file for the compiled routine
error('spm_add.c not compiled - see Makefile')
