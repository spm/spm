function mesh = spm_eeg_inv_transform_mesh(M, mesh)
% Applies homogenous transformation to cortical mesh
% FORMAT mesh = spm_eeg_inv_transform_mesh(M, mesh)
%
% M1          - affine transformation matrix [4 x 4]
% mesh        - 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_mesh.m 2863 2009-03-11 20:25:33Z guillaume $

fn = fieldnames(mesh);

tess_ind = strmatch('tess_', fn);
tess_ind = setdiff(tess_ind, strmatch('tess_mni', fn, 'exact'));

for i = 1:length(tess_ind)
    cmesh =  export(gifti(getfield(mesh, fn{tess_ind(i)})), 'spm');
    cmesh.vert = spm_eeg_inv_transform_points(M,cmesh.vert);
    mesh = setfield(mesh, fn{tess_ind(i)}, cmesh);
end
