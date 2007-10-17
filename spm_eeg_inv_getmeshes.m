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
% Modified by Rik if want only, eg, skull, scalp meshes
% $Id: spm_eeg_inv_getmeshes.m 954 2007-10-17 15:12:26Z rik $

if nargin == 2
  num_mesh = varargin{2};
  if num_mesh<1 | num_mesh>3
    error('number of meshes must lie between 1 and 3')
  end
else
  num_mesh = 3;		% Do all by default (skull, scalp, cortex)
end

% checks and defaults
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{1});
Msize   = D.inv{val}.mesh.Msize;
Nsize   = [3127 4207 5122 7222];
Nvert   = [2002 2002 Nsize(Msize)];

mesh_labels = strvcat('Tess. inner skull','Tess. scalp','Tess. cortex');
Il{1}   = D.inv{val}.mesh.msk_iskull;
Il{2}   = D.inv{val}.mesh.msk_scalp;
Il{3}   = D.inv{val}.mesh.msk_cortex;

[Centre_vx,Centre_mm]     = spm_eeg_inv_CtrBin(Il{1});
D.inv{val}.mesh.Centre_vx = Centre_vx;
D.inv{val}.mesh.Centre_mm = Centre_mm;

% create meshes
%==========================================================================
fprintf('%c','='*ones(1,80)), fprintf('\n')

for m = 1:num_mesh
 fprintf(['\tGenerate %s mesh from binary volume\n'],mesh_labels(m,:));
 head(m) = spm_eeg_inv_TesBin(Nvert(m),Centre_mm,Il{m},mesh_labels(m,:));
end

% smooth
%--------------------------------------------------------------------------
for m = 1:num_mesh
 fprintf(['\tRegularising %s vertices\n'],mesh_labels(m,:));
 head(m) = spm_eeg_inv_ElastM(head(m));
end

% skull mesh
%--------------------------------------------------------------------------
vert = head(1).XYZmm';
face = head(1).tri';
norm = spm_eeg_inv_normals(vert,face);
D.inv{val}.mesh.tess_iskull.vert = vert;
D.inv{val}.mesh.tess_iskull.face = face;
D.inv{val}.mesh.tess_iskull.norm = norm;
D.inv{val}.mesh.Iskull_Nv        = length(vert);
D.inv{val}.mesh.Iskull_Nf        = length(face);

if num_mesh > 1
 % scalp mesh
 %--------------------------------------------------------------------------
 vert = head(2).XYZmm';
 face = head(2).tri';
 norm = spm_eeg_inv_normals(vert,face);
 D.inv{val}.mesh.tess_scalp.vert = vert;
 D.inv{val}.mesh.tess_scalp.face = face;
 D.inv{val}.mesh.tess_scalp.norm = norm;
 D.inv{val}.mesh.Scalp_Nv        = length(vert);
 D.inv{val}.mesh.Scalp_Nf        = length(face);

% cortex mesh
%--------------------------------------------------------------------------
 if num_mesh > 2
  vert = head(3).XYZmm';
  face = head(3).tri';
  norm = spm_eeg_inv_normals(vert,face);
  D.inv{val}.mesh.tess_ctx.vert  = vert;
  D.inv{val}.mesh.tess_ctx.face  = face;
  D.inv{val}.mesh.tess_ctx.norm  = norm;
  D.inv{val}.mesh.Ctx_Nv         = length(vert);
  D.inv{val}.mesh.Ctx_Nf         = length(face);
 end
end

varargout{1} = D;
fprintf('%c','='*ones(1,80)), fprintf('\n')
return
