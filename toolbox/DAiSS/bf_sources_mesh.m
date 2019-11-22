function mesh = bf_sources_mesh(BF, S)
% Generate cortical mesh
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources_mesh.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    orient         = cfg_menu;
    orient.tag     = 'orient';
    orient.name    = 'How to orient the sources';
    orient.labels  = {'Unoriented', 'Original', 'Downsampled'};
    orient.values  = {'unoriented', 'original', 'downsampled'};
    orient.val     = {'unoriented'};
    
    fdownsample      = cfg_entry;
    fdownsample.tag  = 'fdownsample';
    fdownsample.name = 'Downsample factor';
    fdownsample.strtype = 'r';
    fdownsample.num = [1 1];
    fdownsample.val = {1};
    fdownsample.help = {'A number that determines mesh downsampling',...
        'e.g 5 for taking every 5th vertex'};
    
    symmetric         = cfg_menu;
    symmetric.tag     = 'symmetric';
    symmetric.name    = 'Symmetric';
    symmetric.help    = {'Create a symmetric mesh by reflecting on of the hemispheres.'};
    symmetric.labels  = {'No', 'Reflect left', 'Reflect right'};
    symmetric.values  = {'no', 'left', 'right'};
    symmetric.val = {'no'};
    
    flip         = cfg_menu;
    flip.tag     = 'flip';
    flip.name    = 'Flip';
    flip.help    = {'Flip the mesh relative to midsagittal plane.'};
    flip.labels  = {'yes', 'no'};
    flip.values  = {true, false};
    flip.val = {false};
    
    mesh = cfg_branch;
    mesh.tag = 'mesh';
    mesh.name = 'Cortical mesh';
    mesh.val = {orient, fdownsample, symmetric, flip};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

original = BF.data.mesh.tess_mni;

canonical = original;

if S.fdownsample ~= 1
    canonical = export(gifti(reducepatch(export(gifti(canonical), 'patch'), 1/S.fdownsample)), 'spm');
end

if ~isequal(S.symmetric, 'no')
    meshcomp = spm_mesh_split(gifti(canonical));
    
    if isequal(S.symmetric, 'left')
        if median(meshcomp(1).vertices(:, 1))<0
            side = meshcomp(1);
        else
            side = meshcomp(2);
        end
    elseif isequal(S.symmetric, 'right')
        if median(meshcomp(1).vertices(:, 1))>0
            side = meshcomp(1);
        else
            side = meshcomp(2);
        end
    end
    
    cside = side;
    cside.vertices(:, 1) = -1*side.vertices(:, 1);
    cside.faces = side.faces+size(side.vertices, 1);
    
    canonical = [];
    if isequal(S.symmetric, 'left')
        canonical.vert = [side.vertices; cside.vertices];
    else
        canonical.vert = [cside.vertices; side.vertices];
    end
    
    canonical.face = [side.faces; cside.faces];
end

mesh = [];
mesh.canonical = canonical;


if isfield(BF.data.mesh, 'def')
    mesh.individual = spm_swarp(gifti(mesh.canonical), BF.data.mesh.def);
    original        = spm_swarp(gifti(original), BF.data.mesh.def);
else
    mesh.individual = mesh.canonical;
end

M = BF.data.transforms.toNative;

mesh.individual      = export(gifti(mesh.individual), 'spm');
mesh.individual.vert = spm_eeg_inv_transform_points(inv(M), mesh.individual.vert);

original = export(gifti(original), 'spm');
original.vert = spm_eeg_inv_transform_points(inv(M), original.vert);

mesh.pos = mesh.individual.vert;

if S.flip
    M1 = eye(4);
    M1(1, 1) = -1;
    M1 = BF.data.transforms.toMNI_aligned\M1*BF.data.transforms.toMNI_aligned;
    mesh.pos = spm_eeg_inv_transform_points(M1, mesh.individual.vert);
end

switch S.orient
    case 'original'
        norm = spm_mesh_normals(export(gifti(original), 'patch'), true);
        if S.fdownsample == 1
            mesh.ori = norm;
        else
            mesh.ori = 0*mesh.pos;
            for i = 1:size(mesh.pos, 1)
                mesh.ori(i, :) = norm(all(repmat(mesh.pos(i, :), size(original.vert, 1), 1) == original.vert, 2), :);
            end
        end
    case 'downsampled'
        mesh.ori = spm_mesh_normals(export(gifti(mesh.individual), 'patch'), true);
end
