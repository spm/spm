function GX = spm_global(V)
% returns the global mean for a volume image - a compiled routine
% FORMAT GX = spm_global(V)
% V   - memory mapped volume
% GX  - mean global activity
%__________________________________________________________________________
%
% spm_global returns the mean counts integrated over all the  
% slices from the volume.
%
% The mean is estimated after discounting voxels outside the object
% using a criteria of greater than > (global mean)/8.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Anonymous
% $Id: spm_global.m 4182 2011-02-01 12:29:09Z guillaume $


%-This is merely the help file for the compiled routine
error('spm_global.c not compiled - see Makefile')

