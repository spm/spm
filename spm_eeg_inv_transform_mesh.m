function mesh = spm_eeg_inv_transform_mesh(M, mesh)
% Applies affine transformation to surface mesh
% FORMAT mesh = spm_eeg_inv_transform_mesh(M, mesh)
%
% M           - affine transformation matrix [4 x 4]
% mesh        - patch structure
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_mesh.m 4082 2010-10-07 16:11:48Z guillaume $

fn = fieldnames(mesh);

tess_ind = find(strncmp(fn,'tess_',5) & ~strcmp(fn,'tess_mni'));

for i = 1:length(tess_ind)
    cmesh = export(gifti(mesh.(fn{tess_ind(i)})),'spm');
    cmesh.vert = spm_eeg_inv_transform_points(M,cmesh.vert);
    mesh.(fn{tess_ind(i)}) = cmesh;
end
