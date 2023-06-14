function F = spm_mesh_sdf(M, V, m)
% Compute the signed distance field (SDF) to a triangle mesh
% FORMAT D = spm_mesh_sdf(M, V, m)
% M        - a patch structure with fields 'faces' and 'vertices'
% V        - an spm_vol structure with fields 'dim' and 'mat'
% m        - a binary mask (image filename or spm_vol structure)
%            [default: none]
%
% F        - a 3D array containing signed distance values
%__________________________________________________________________________
%
% Example:
%
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
% M = export(M,'patch');
% M.faces = double(M.faces);
% V = spm_vol(fullfile(spm('Dir'),'canonical','single_subj_T1.nii'));
%
% F = spm_mesh_sdf(M, V);
%
% D = struct(...
%     'fname',   'sdf.nii',...
%     'dim',     V.dim,...
%     'dt',      [spm_type('float64') spm_platform('bigend')],...
%     'mat',     V.mat,...
%     'pinfo',   [1 0 0]');
% spm_write_vol(D, F);
%
% spm_check_registration(D.fname);
% spm_ov_mesh('Display', 1, M);
% spm_colourmap('hot')
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if nargin < 3
    m = [];
else
    m = spm_vol(m);
    spm_check_orientations(...
        struct('dim',{V.dim,m.dim}, 'mat',{V.mat,m.mat}));
    m = spm_read_vols(m);
    m = reshape(logical(m),1,[]);
end

I = 1:prod(V.dim);
if ~isempty(m)
    I(~m) = [];
end
[X,Y,Z] = ind2sub(V.dim, I');
XYZ = spm_mesh_transform([X,Y,Z], V.mat);

F = spm_mesh_dist(M, XYZ);

if isempty(m)
    F = reshape(F, V.dim);
else
    f = NaN(V.dim);
    f(m) = F;
    F = f;
end
