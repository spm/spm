function [vol, fid, mesh] = spm_eeg_inv_meshing(sMRI, Msize, modality)
% Apply the inverse spatial deformation to the template mesh
% to obtain the individual cortical mesh
% save the individual .mat tesselation of the chosen size
%
% FORMAT [megvol, fid, mesh] = spm_eeg_inv_meshing(filename, Msize)
% Input:
% sMRI - name of the sMRI file
% Msize - size of the mesh (1-4)
% Output:
% megvol - single shell BEM for MEG
% fid    - fiducials (head surface + points inverse normalized from the template)
% mesh   - inverse - normalized canonical mesh
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_meshing.m 1726 2008-05-26 16:45:55Z vladimir $

if nargin <2
    Msize = 4;
end

if nargin <3
    modality = 'MEG';
end

%%
mesh = [];
mesh.Msize = Msize;
mesh.template = 0;
mesh.sMRI = sMRI;

% Segment and normalise structural
%======================================================================
[mesh]  = spm_eeg_inv_segment(mesh);

% Compute the masks
%======================================================================
mesh  = spm_eeg_inv_getmasks(mesh);
%%
% Get volume models and head surface meshes
%------------------------------------------------------------------
[mesh, vol, fid] = spm_eeg_inv_getmeshes(mesh, modality);
%%
% Canonical cortical mesh
%------------------------------------------------------------------
[junk, template_fid, template_mesh] = spm_eeg_inv_template(Msize, modality);

mesh.tess_mni.vert    = template_mesh.tess_mni.vert;
mesh.tess_mni.face    = template_mesh.tess_mni.face; 

mesh.tess_ctx.vert    = mesh.tess_mni.vert;
mesh.tess_ctx.face    = mesh.tess_mni.face ;

% Inverse transform the template mesh
%------------------------------------------------------------------
mesh = spm_eeg_inv_transform_mesh(inv(mesh.Affine), mesh);

% Get labeled fiducial points by inverse transforming template fiducials
%------------------------------------------------------------------
fid.fid = getfield(forwinv_transform_headshape(inv(mesh.Affine), template_fid), 'fid');

spm('Pointer','Arrow');
