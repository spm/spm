function varargout = spm_get_lm(varargin)
% Identification of local maxima in 3(or 2)D volume - a compiled routine
%
% FORMAT: INDEX = spm_get_lm(VOL,LIST)
%
% Routine that identifies which voxels in a list of coordinates
% that are local maxima, and returns a list of indicies into
% the coordinate list for those maxima.
%
% Input:
% VOL          : 3(or 2)D volume of statistics (e.g. t or F)
% LIST         : 3xn (or 2xn) list of voxel coordinates of 
%                tentative local maxima.
%
% Output:
% INDEX        : Index into LIST such that LIST(:,INDEX)
%                returns those coordinates that are truly
%                local maxima.
%_______________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: spm_get_lm.m 2696 2009-02-05 20:29:48Z guillaume $

%-This is merely the help file for the compiled routine
error('spm_get_lm.c not compiled - see Makefile');
