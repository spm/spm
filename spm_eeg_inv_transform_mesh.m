function mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% Applies homogenous transformation to cortical mesh.
% Loads the mesh from gifti files if necessary
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_mesh.m 2720 2009-02-09 19:50:46Z vladimir $

fn = fieldnames(mesh);

tess_ind = strmatch('tess_', fn);
tess_ind = setdiff(tess_ind, strmatch('tess_mni', fn, 'exact'));

for i = 1:length(tess_ind)
    cmesh =  export(gifti(getfield(mesh, fn{tess_ind(i)})), 'spm');
    old = cmesh.vert;
    old(:,4) = 1;
    new = old * M1';
    new = new(:,1:3);
    cmesh.vert = new;
    mesh = setfield(mesh, fn{tess_ind(i)}, cmesh);
end
