function mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% mesh = spm_eeg_inv_transform_mesh(M1, mesh)
% Applies homogenous transformation to cortical mesh. 
% _________________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_transform_mesh.m 1523 2008-04-30 17:33:04Z vladimir $

old = mesh.tess_ctx.vert;
old(:,4) = 1;
new = old * M1';
new = new(:,1:3);
mesh.tess_ctx.vert = new;