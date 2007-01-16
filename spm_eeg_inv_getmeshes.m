function varargout = spm_eeg_inv_getmeshes(varargin);

%==========================================================================
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_getmeshes(D)
% Input:
% D		    - input data struct (optional)
% Output:
% D			- data struct with the tessellation fields
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_getmeshes.m 716 2007-01-16 21:13:50Z karl $


% checks and defaults
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
Iisk    = D.inv{val}.mesh.msk_iskull;
Ictx    = D.inv{val}.mesh.msk_iskull;
Iscl    = D.inv{val}.mesh.msk_scalp;
Msize   = D.inv{val}.mesh.Msize;
Nsize   = [3127 4207 5122 7222]
Nvert   = [2002 2002 Nsize(Msize)];

[Centre_vx,Centre_mm]     = spm_eeg_inv_CtrBin(Iisk);
D.inv{val}.mesh.Centre_vx = Centre_vx;
D.inv{val}.mesh.Centre_mm = Centre_mm;

% create meshes
%==========================================================================
fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
fprintf(['\tGenerate surface meshes from binary volumes\n']);

mesh_labels = strvcat('Tess. inner skull','Tess. scalp','Tess. cortex');
head(1) = spm_eeg_inv_TesBin(Nvert(1),Centre_mm,Iisk,mesh_labels(1,:));
head(2) = spm_eeg_inv_TesBin(Nvert(2),Centre_mm,Iscl,mesh_labels(2,:));
head(3) = spm_eeg_inv_TesBin(Nvert(3),Centre_mm,Ictx,mesh_labels(3,:));

% smooth
%--------------------------------------------------------------------------
fprintf(['\t\tRegularising vertices\n']);
head(1) = spm_eeg_inv_ElastM(head(1));
head(2) = spm_eeg_inv_ElastM(head(2));
head(3) = spm_eeg_inv_ElastM(head(3));

% skull mesh
%--------------------------------------------------------------------------
vert = head(1).XYZmm';
face = head(1).tri';
norm = spm_eeg_inv_normals(vert,face);
D.inv{val}.mesh.tess_iskull.vert = vert;
D.inv{val}.mesh.tess_iskull.vert = face;
D.inv{val}.mesh.tess_iskull.vert = norm;
D.inv{val}.mesh.Iskull_Nv        = length(vert);
D.inv{val}.mesh.Iskull_Nf        = length(face);

% scalp mesh
%--------------------------------------------------------------------------
vert = head(2).XYZmm';
face = head(2).tri';
norm = spm_eeg_inv_normals(vert,face);
D.inv{val}.mesh.tess_scalp.vert = vert;
D.inv{val}.mesh.tess_scalp.vert = face;
D.inv{val}.mesh.tess_scalp.vert = norm;
D.inv{val}.mesh.Scalp_Nv        = length(vert);
D.inv{val}.mesh.Scalp_Nf        = length(face);

% cortex mesh
%--------------------------------------------------------------------------
vert = head(3).XYZmm';
face = head(3).tri';
norm = spm_eeg_inv_normals(vert,face);
D.inv{val}.mesh.tess_ctx.vert  = vert;
D.inv{val}.mesh.tess_ctx.vert  = face;
D.inv{val}.mesh.tess_ctx.vert  = norm;
D.inv{val}.mesh.Ctx_Nv         = length(vert);
D.inv{val}.mesh.Ctx_Nf         = length(face);

varargout{1} = D;
fprintf('%c','='*ones(1,80)), fprintf('\n')
return
