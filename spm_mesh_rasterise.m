function R = spm_mesh_rasterise(M, V)
% Rasterise a triangle mesh on a regular grid
% FORMAT R = spm_mesh_rasterise(M, V)
% M        - a patch structure or GIfTI object
% V        - structure with fields 'dim' and 'mat' defining the grid
%
% R        - logical array: 1 for inside and 0 for outside
%__________________________________________________________________________
%
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii'));
% V = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
% R = spm_mesh_rasterise(M, V);
% V.fname = 'raster.nii';
% V.dt(1) = spm_type('uint8');
% V.pinfo = [1 0 0]';
% V.dat = uint8(R);
% spm_check_registration(V)
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


I = 1:prod(V.dim);
[X,Y,Z] = ind2sub(V.dim, I');
XYZ = spm_mesh_transform([X,Y,Z], V.mat);

R = false(V.dim);
for i=1:size(XYZ,1)
    ray = struct('orig',XYZ(i,:)', 'vec',XYZ(i,:)' + [0 0 1]');
    R(i) = mod(nnz(spm_mesh_ray_intersect(M, ray)), 2);
end
