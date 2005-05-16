function varargout = spm_get_lm(varargin)
%
% Identification of local maxima in 3(or 2)D volume. 
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jesper Andersson 
% $Id: spm_get_lm.m 159 2005-05-16 14:00:56Z guillaume $


error('spm_get_lm.c not compiled - see Makefile');
