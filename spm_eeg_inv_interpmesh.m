function H = spm_eeg_inv_interpmesh(S)

%=======================================================================
% Compute the interpolation matrix which enables to infer voxel values from
% from vertices values on the vertices.
% This computation is based on a trilinear interpolation which is done by
% calling the BrainSTorm 2.0 bst_tri_interp.m function
% ( see http://neuroimage.usc.edi/brainstorm )
%
% FORMAT H = spm_eeg_inv_interpmesh(S)
% Input:
% S		    - input data struct (optional)
% Output:
% H			- interpolation matrix from 3D mesh to voxel space
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_interpmesh.m 621 2006-09-12 17:22:42Z karl $

spm_defaults

try
    D = S;
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D = spm_eeg_ldata(D);
end

try
    val = D.val;
catch
    val = length(D.inv);
end

% Load the mesh
load(D.inv{val}.mesh.tess_ctx);

% Converting Vertices coordinates into indices in the MRI volume
Vo      = spm_vol(D.inv{val}.mesh.sMRI);
dims    = Vo.dim(1:3);
pixdim  = Vo.private.hdr.dime.pixdim([2 3 4]);
vert    = round(vert(:,1:3)./(ones(length(vert),1)*pixdim));

% Interpolating the Mesh into the 3D MRI voxel space (BST function)
H       = bst_tri_interp(dims,vert,face);
