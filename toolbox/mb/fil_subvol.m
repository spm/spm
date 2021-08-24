function [d,M] = fil_subvol(Nii,bb)
% Dimensions and voxel-to world mapping of a subvolume
% FORMAT VO = subvol(Nii,bb)
% V      - SPM NIfTI object
% bb     - bounding box (2 x 3)
%__________________________________________________________________________
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% John Ashburner
% $Id: fil_subvol.m 8139 2021-08-24 19:38:01Z guillaume $

if ~isa(Nii,'nifti')
    Nii = nifti(Nii);
end
bb0 =  [1 1 1; size(Nii.dat,1) size(Nii.dat,2) size(Nii.dat,3)];
if nargin<2
    bb = bb0;
else
    if ~all(size(bb)==[2 3])
        error('Incorrectly sized bb.');
    end
    msk     = ~isfinite(bb);
    bb(msk) = bb0(msk);
    bb      = sort(bb);
    bb(1,:) = round(max(bb(1,:),bb0(1,:)));
    bb(2,:) = round(min(bb(2,:),bb0(2,:)));
end
d = diff(bb)+1;
M = Nii.mat*spm_matrix((bb(1,:)-1));
%==========================================================================

