function varargout = spm_map_vol(varargin)
% memory map of a volume image - compiled routine
% FORMAT V = spm_map_vol(P,[DIM VOX SCALE TYPE OFFSET]);
% P      - filename
% DIM    -  [x y z] - image size {voxels}
% VOX    -  [x y z] - voxel size {mm}
% SCALE  -  scaling coefficient to multiply voxel values by
% TYPE   -  data type  (see spm_type.m for supported types and specifiers)
% OFFSET -  offset of the first byte of the volume (bytes)
%_______________________________________________________________________
%
% spm_map_vol returns a vector V identifying a memory mapped image
% volumne on disk.  Memory mapping avoids having very large objects
% in working memory
%
% see also spm_map.m and spm_unmap_vol.m
%
%_______________________________________________________________________
% %W% John Ashburner %E%

%-This is merely the help file for the compiled routine
error('spm_map_vol.c not compiled - see spm_MAKE.sh')
