function F = spm_mesh_sdf(M, V)
% Compute the signed distance field (SDF) to a triangle mesh
% FORMAT D = spm_mesh_sdf(M, V)
% M        - a patch structure with fields 'faces' and 'vertices'
% V        - an spm_vol structure with fields 'dim' and 'mat'
%
% F        - a 3D array containing signed distance values
%__________________________________________________________________________
%
% Example:
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
% spm_write_vol(O,D);
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


I = 1:prod(V.dim);
[X,Y,Z] = ind2sub(V.dim, I');
XYZ = spm_mesh_transform([X,Y,Z], V.mat);

F = spm_mesh_dist(M, XYZ);
F = reshape(F, V.dim);
