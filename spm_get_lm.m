function varargout = spm_get_lm(varargin)
%
% Identification of local maxima in 3D volume. 
%
% FORMAT: INDEX = spm_bwlabel(VOL,LIST)
%
% Routine that identifies which voxels in a list of coordinates
% that are local maxima, and returns a list of indicies into
% the coordinate list for those maxima.
%
% Input:
% VOL          : 3D volume of statistics (e.g. t or F)
% LIST         : 3xn list of voxel coordinates of 
%                tentative local maxima.
%
% Output:
% INDEX        : Index into LIST such that LIST(:,INDEX)
%                returns those coordinates that are truly
%                local maxima.
%_______________________________________________________________
% %W% Jesper Andersson %E% 

error('spm_get_lm.c has not been compiled');
