
% compute residual variance - a compiled routine
% FORMAT s = spm_resid2(VI,VO,U)
% VI    - Vector of mapped volumes (from spm_map or spm_vol).
% VO    - Description of output volume that gets passed to
%         spm_write_plane.m
% U     - orthogonalized design matrix.
% s     - Scalefactor for output image.
%_______________________________________________________________________
%
% Each voxel of the output volume is computed by sum((x-U*(U'*x)).^2),
% where U is the orthogonalized disign matrix, and x is a time series
% of voxel intensities.
% The formula for computing the residuals obtained using a design matrix
% A is normally sum((x-A*(A\x)).^2).  Decomposing A using singular value
% decomposition (SVD) gives A = U*S*V' and A^(-1) = V'*(S^(-1))*U.
% Therefore, x-A*(A\x) = x-U*S*V'*(V'*(S^(-1))*U*x) = x-U*(U'*x).
%
% See also, spm_add.m.
%_______________________________________________________________________
% %W% John Ashburner %E%
