function [mesh, vol, fid] = spm_eeg_inv_getmeshes(mesh, modality)
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT D = spm_eeg_inv_getmeshes(D)
% Input:
%   mesh      - input data struct with the binary mask volumes filenames
%   modality  - the kind of BEM to make - EEG or MEG
% Output:
%   mesh      - mesh struct with the tessellation fields
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout, Christophe Phillips, Rik Henson
% $Id: spm_eeg_inv_getmeshes.m 1727 2008-05-26 17:49:22Z vladimir $


if nargin<2
    modality = 'MEG';
end

switch modality
    case 'MEG'
        lbuild = [1,3];
    case 'EEG'
        lbuild =  [1:3];
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

% prepare forwinv volume models from meshes
%==========================================================================
% bnd = [];
% scalp mesh
%--------------------------------------------------------------------------
switch modality
    case 'MEG'
        ind = find(list==3);
        bnd = struct('pnt',head(ind).XYZmm', 'tri',head(ind).tri');
        vol = struct('bnd', bnd, 'type', 'nolte');
    case 'EEG'
        for i = 1:3
            ind = find(list==i);
            bnd(i) = struct('pnt',head(ind).XYZmm', 'tri',head(ind).tri');
        end
        vol = struct('bnd',bnd,'cond',[0.3300 0.0041 0.3300],'type','dipoli');
end

fid = [];
fid.pnt = head(find(list==1)).XYZmm';
fid.tri = head(find(list==1)).tri';
fid.fid.pnt = [];
fid.fid.label = {};

vol = forwinv_convert_units(vol, 'mm');
fid = forwinv_convert_units(fid, 'mm');


% cortex mesh
%--------------------------------------------------------------------------
% if any(list==4)
%     ind = find(list==4);
%     vert = head(ind).XYZmm';
%     face = head(ind).tri';
%     norm = spm_eeg_inv_normals(vert,face);
%     mesh.tess_ctx.vert  = vert;
%     mesh.tess_ctx.face  = face;
%     mesh.tess_ctx.norm  = norm;
%     mesh.Ctx_Nv         = length(vert);
%     mesh.Ctx_Nf         = length(face);
% end

fprintf('%c','='*ones(1,80)), fprintf('\n')
return
