
% averages a series of images - a compiled routine
% FORMAT s = spm_mean(n,TYPE,Q,P,S)
% n	-	number of voxels per image
% TYPE	- 	data type (see spm_type.m & volume.h)
% Q	-	filename for averaged image
% P	-	matrix of filenames to be averaged (rowwise strings)
% S	-	optional vector of scalefactors / weights
% s	-	the scalefactor which should be assigned to Q
%_______________________________________________________________________
%
% spm_mean computes a weighted sum of a set of image files to produce
% an mean image that is written to a named file (Q). No headers are
% read or written, in particular scalefactors are ignored. The image is
% written in the same type as the input images.
%
% The weights (S) default to 1/size(P,1) - resulting an average of the
% image files being written to Q. For a "proper" average including
% scalefactors, the weights s should be specified as sf/length(sf) for
% sf a vector of scalefactors for the image files specified in P.
%
% See also: spm_average - for details on how to usefully use this function
%
%_______________________________________________________________________
% %W% Jean-Baptiste Poline, John Ashburner%E%
