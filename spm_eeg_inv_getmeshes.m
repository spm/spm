function [mesh,vol] = spm_eeg_inv_getmeshes(mesh,lbuild)
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_getmeshes(D)
% Input:
%   mesh      - input data struct with the binary mask volumes filenames
%   lbuild    - list of meshes to build (scalp 1, oskull 2, iskull 3, cortex 4)
% Output:
%   mesh      - mesh struct with the tessellation fields
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Christophe Phillips, Rik Henson
% $Id: spm_eeg_inv_getmeshes.m 1488 2008-04-27 14:11:48Z vladimir $

if nargin<1
    D = spm_eeg_load;
    [D,val] = spm_eeg_inv_check(varargin{1});
    mesh = D.inv{val}.mesh;
end
if nargin<2
    lbuild = 1:4;
end

% checks and defaults
%--------------------------------------------------------------------------
list = 1:4;
Msize = mesh.Msize;
Nsize = [3127 4207 5122 7222];
Nvert   = [502 1622 2002 Nsize(Msize)];
% By default 502/1622/2002 vertices per "volume surface" (scalp/oskull/iskull)
% for the cortex, use the "usual" nr of vertices.
% NOTE:
% I wouldn't trust too much the tessalation of the cortical surface, as
% we're NOT using any sophisticated surface extraction procedure like in
% BrainVisa/BrainVoyager...

% Check what should be built
%--------------------------------------------------------------------------
Il = {};
mesh_labels = strvcat('Tess. scalp','Tess. outer skull', ...
                        'Tess. inner skull','Tess. cortex');
if ~isempty(mesh.msk_scalp) && any(lbuild==1)
    Il{end+1}   = mesh.msk_scalp;
else
    ind = find(list==1);
    mesh_labels(ind,:) = [];
    list(ind) = [];
    Nvert(ind) = [];
end
if ~isempty(mesh.msk_oskull) && any(lbuild==2)
    Il{end+1}   = mesh.msk_oskull;
else
    ind = find(list==2);
    mesh_labels(ind,:) = [];
    list(ind) = [];
    Nvert(ind) = [];
end
if ~isempty(mesh.msk_iskull) && any(lbuild==3)
    Il{end+1}   = mesh.msk_iskull;
else
    ind = find(list==3);
    mesh_labels(ind,:) = [];
    list(ind) = [];
    Nvert(ind) = [];
end
if ~isempty(mesh.msk_cortex) && any(lbuild==4)
    Il{end+1}   = mesh.msk_cortex;
else
    ind = find(list==4);
    mesh_labels(ind,:) = [];
    list(ind) = [];
    Nvert(ind) = [];
end

if isempty(list)
    error('No binary mask available ! Can''t create meshes...')
end

[mesh.Centre_vx,mesh.Centre_mm]      = spm_eeg_inv_CtrBin(Il{end});

% create meshes
%==========================================================================
fprintf('%c','='*ones(1,80)), fprintf('\n')
num_mesh = length(list);
for m = 1:num_mesh
    fprintf(['\tGenerate %s mesh from binary volume\n'],mesh_labels(m,:));
    head(m) = spm_eeg_inv_TesBin(Nvert(m),mesh.Centre_mm,Il{m},mesh_labels(m,:));
end

% smooth
%--------------------------------------------------------------------------
for m = 1:num_mesh
    fprintf(['\tRegularising %s vertices\n'],mesh_labels(m,:));
    head(m) = spm_eeg_inv_ElastM(head(m));
end

% store meshes
%==========================================================================
bnd = [];
% scalp mesh
%--------------------------------------------------------------------------
if any(list==1)
    ind = find(list==1);
    bnd(end+1) = struct('pnt',head(ind).XYZmm', ...
                        'tri',head(ind).tri');
%     vert = head(ind).XYZmm';
%     face = head(ind).tri';
%     norm = spm_eeg_inv_normals(vert,face);
%     mesh.tess_scalp.vert  = vert;
%     mesh.tess_scalp.face  = face;
%     mesh.tess_scalp.norm  = norm;
%     mesh.Iskull_Nv        = length(vert);
%     mesh.Iskull_Nf        = length(face);
end

% outer skull mesh
%--------------------------------------------------------------------------
if any(list==2)
    ind = find(list==2);
    bnd(end+1) = struct('pnt',head(ind).XYZmm', ...
                        'tri',head(ind).tri');
%     vert = head(ind).XYZmm';
%     face = head(ind).tri';
%     norm = spm_eeg_inv_normals(vert,face);
%     mesh.tess_oskull.vert = vert;
%     mesh.tess_oskull.face = face;
%     mesh.tess_oskull.norm = norm;
%     mesh.Oskull_Nv        = length(vert);
%     mesh.Oskull_Nf        = length(face);
end

% inner skull mesh
%--------------------------------------------------------------------------
if any(list==3)
    ind = find(list==3);
    bnd(end+1) = struct('pnt',head(ind).XYZmm', ...
                        'tri',head(ind).tri');
%     vert = head(ind).XYZmm';
%     face = head(ind).tri';
%     norm = spm_eeg_inv_normals(vert,face);
%     mesh.tess_iskull.vert = vert;
%     mesh.tess_iskull.face = face;
%     mesh.tess_iskull.norm = norm;
%     mesh.Scalp_Nv         = length(vert);
%     mesh.Scalp_Nf         = length(face);
end

vol = struct('bnd',bnd,'cond',[0.3300 0.0041 0.3300],'type','dipoli');

% cortex mesh
%--------------------------------------------------------------------------
if any(list==4)
    ind = find(list==4);
    vert = head(ind).XYZmm';
    face = head(ind).tri';
    norm = spm_eeg_inv_normals(vert,face);
    mesh.tess_ctx.vert  = vert;
    mesh.tess_ctx.face  = face;
    mesh.tess_ctx.norm  = norm;
    mesh.Ctx_Nv         = length(vert);
    mesh.Ctx_Nf         = length(face);
end

fprintf('%c','='*ones(1,80)), fprintf('\n')
return
