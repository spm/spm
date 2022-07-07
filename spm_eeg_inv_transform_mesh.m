function mesh = spm_eeg_inv_transform_mesh(M, mesh)
% Apply affine transformation to surface mesh
% FORMAT mesh = spm_eeg_inv_transform_mesh(M, mesh)
%
% M           - affine transformation matrix [4 x 4]
% mesh        - patch structure
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


fn = fieldnames(mesh);

tess_ind = find(strncmp(fn,'tess_',5) & ~strcmp(fn,'tess_mni'));

for i = 1:length(tess_ind)
    cmesh = export(gifti(mesh.(fn{tess_ind(i)})),'spm');
    cmesh.vert = spm_eeg_inv_transform_points(M,cmesh.vert);
    mesh.(fn{tess_ind(i)}) = cmesh;
end

if isfield(mesh, 'fid')
    mesh.fid.pnt = spm_eeg_inv_transform_points(M, mesh.fid.pnt);
    mesh.fid.fid.pnt = spm_eeg_inv_transform_points(M, mesh.fid.fid.pnt);
end
