function [d,M] = fil_subvol(Nii,bb)
% Dimensions and voxel-to world mapping of a subvolume
% FORMAT [d,M] = fil_subvol(Nii,bb)
% Nii    - SPM NIfTI object
% bb     - bounding box (2 x 3)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2021-2022 Wellcome Centre for Human Neuroimaging


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
