function mesh = bf_sources_scalp(BF, S)
% Generate source space on the scalp surface, as a part of measuring a
% magnetomyogram (MMG) of the face.
% 
% See https://doi.org/10.1111/psyp.13507 for more information.
%__________________________________________________________________________


% George O'Neill
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

if nargin == 0
    
    orient         = cfg_menu;
    orient.tag     = 'orient';
    orient.name    = 'How to orient the sources';
    orient.labels  = {'Unoriented', 'Original'};
    orient.values  = {'unoriented', 'original'};
    orient.val     = {'unoriented'};
    
%     warp           = cfg_menu;
%     warp.tag       = 'warp';
%     warp.name      = 'How to warp between canonical and individual mesh';
%     warp.labels    = {'Affine Transform', 'Deformation Field'};
%     warp.values    = {'affine','def'};
%     warp.val       = {'affine'};
%     
    mesh = cfg_branch;
    mesh.tag = 'scalp';
    mesh.name = 'Scalp mesh';
    mesh.val = {orient};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

% if ~isfield(S,'warp'), S.warp = 'affine'; end

original = BF.data.mesh.tess_scalp;
istemplate = BF.data.mesh.template;

mesh.individual = original;

M = BF.data.transforms.toNative;
A = BF.data.mesh.Affine;

mesh.individual      = export(gifti(mesh.individual), 'spm');
mesh.canonical       = mesh.individual;
mesh.individual.vert = spm_eeg_inv_transform_points(inv(M), mesh.individual.vert);
mesh.canonical.vert = spm_eeg_inv_transform_points(A, mesh.canonical.vert);

original = export(gifti(original), 'spm');
original.vert = spm_eeg_inv_transform_points(inv(M), original.vert);

mesh.pos = mesh.individual.vert;

switch S.orient
    case 'original'
        norm = spm_mesh_normals(export(gifti(original), 'patch'), true);
        mesh.ori = norm;
        
end
