function V = spm_imwrite(varargin)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_imwrite(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
% Y          - a one, two or three dimensional matrix containing the image voxels
% V (output) - data structure after modification for writing.
%_______________________________________________________________________
%
% spm_imwrite.m was renamed spm_write_vol.m, to mimic the existing naming
% conventions for memory mapped file handling in SPM.             -andrew
%
%_______________________________________________________________________
% %W% John Ashburner %E%

%-Print warning of obsolescence
%-----------------------------------------------------------------------
warning('spm_imwrite has been renamed spm_write_vo')


%-Pass on arguments to spm_write_vol
%-----------------------------------------------------------------------
V = spm_write_vol(varargin{:})
