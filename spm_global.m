function GX = spm_global(V)
% returns the global mean for a memory mapped volume image
% FORMAT GX = spm_global(V)
% V   - memory mapped volume
% GX  - mean global activity
%_______________________________________________________________________
%
% spm_global returns the mean counts integrated over all the  
% slices from the volume
%
% The mean is estimated after discounting voxels outside the object
% using a criteria of greater than > (global mean)/8
%
%_______________________________________________________________________
% %W% Anon %E%
