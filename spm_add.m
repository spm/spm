
% averages a series of images - a compiled routine
% FORMAT s = spm_mean(n,TYPE,Q,P,S)
% n	-	number of voxels per image
% TYPE	- 	data type - see spm_type.m
% Q	-	filename for averaged image
% P	-	matrix of filenames to be averaged (rowwise strings)
% S	-	optional vector of scalefactors
% s	-	the scalefactor which should be assigned to Q
%____________________________________________________________________________
%
% spm_mean simply averages a set of images to produce an mean image that
% is written to a named file (Q)
%
%__________________________________________________________________________
% %W% %E%
