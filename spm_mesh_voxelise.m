function R = spm_mesh_voxelise(M, V)
% Voxelise a triangle mesh on a regular grid
% FORMAT R = spm_mesh_voxelise(M, V)
% M        - a patch structure or GIfTI object
% V        - structure with fields 'dim' and 'mat' defining the grid
%
% R        - logical array: 1 for inside and 0 for outside
%__________________________________________________________________________
%
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii'));
% V = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
% R = spm_mesh_voxelise(M, V);
% V.fname = 'voxelised.nii';
% V.dt(1) = spm_type('uint8');
% V.pinfo = [1 0 0]';
% V.dat = uint8(R);
% spm_check_registration(V)
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if nargin < 2 || ~isstruct(V)
    if nargin < 2
        vx = 1;
    else
        vx = V;
    end
    if numel(vx) == 1
        vx = [vx vx vx];
    end
    bb = [min(M.vertices,[],1); max(M.vertices,[],1)];
    V =struct('dim', ceil(diff(bb)./vx));
    V.mat = diag([vx 1]);
    V.mat(1:3,4) = bb(1,:) - vx;
end

I = 1:prod(V.dim);
[X,Y,Z] = ind2sub(V.dim, I');
XYZ = spm_mesh_transform([X,Y,Z], V.mat);

R = false(V.dim);

for i=1:size(XYZ,1)
    R(i) = spm_mesh_inside(M,XYZ(i,:));
    fprintf('%d out of %d\n',i,size(XYZ,1));
end
