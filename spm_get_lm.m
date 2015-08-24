function idx = spm_get_lm(vol,list,n)
% Identification of local maxima in 3(or 2)D volume - a compiled routine
% FORMAT idx = spm_get_lm(vol,list,n)
%
% Routine that identifies which voxels in a list of coordinates that are
% local maxima, and returns a list of indices into the coordinate list for
% those maxima.
%
% Input:
% vol      - 3(or 2)D volume of statistics (e.g. t or F)
% list     - 3xn (or 2xn) list of voxel coordinates of tentative local
%            maxima.
% n        - connectivity criterion: 6 (surface), 18 (edge) or 26 (corner).
%            [Default: 18].
%            (for a 2D image these correspond to 4, 8 and 8 respectively).
%
% Output:
% idx      - Index into list such that list(:,idx) returns those
%            coordinates that are truly local maxima.
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: spm_get_lm.m 6534 2015-08-24 16:02:56Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_get_lm.c not compiled - see Makefile');
