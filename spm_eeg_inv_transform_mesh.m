function mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% Applies homogenous transformation to a mesh
% FORMAT mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_mesh.m 2696 2009-02-05 20:29:48Z guillaume $

old = mesh.tess_ctx.vert;
old(:,4) = 1;
new = old * M1';
new = new(:,1:3);
mesh.tess_ctx.vert = new;