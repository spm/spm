function [R, V] = spm_mesh_voxelise(M, V)
% Voxelise a triangle mesh on a regular grid
% FORMAT [R, V] = spm_mesh_voxelise(M, V)
% M        - a patch structure or GIfTI object
% V        - structure with fields 'dim' and 'mat' defining the grid
%            or voxel size for automatic field of view [default: 1]
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


%-Define field of view if only voxel size is provided
%--------------------------------------------------------------------------
if nargin < 2 || ~isstruct(V)
    if nargin < 2
        vx = 1;
    else
        vx = V;
    end
    if numel(vx) == 1
        vx = [vx vx vx];
    end
    bb = spm_mesh_bounding_volume(M);
    V = struct('dim', ceil(diff(bb)./vx));
    V.mat = diag([vx 1]);
    V.mat(1:3,4) = bb(1,:) - vx;
end

%-Compute +z ray intersections from 2D grid at z=1
%--------------------------------------------------------------------------
M = spm_mesh_transform(M, inv(V.mat));

%-Get vertices coordinates for all triangles (for speed)
v1 = M.vertices(M.faces(:,1),:)';
v2 = M.vertices(M.faces(:,2),:)';
v3 = M.vertices(M.faces(:,3),:)';
v  = [v1;v2;v3];

% h = -Inf(V.dim);
% ray = struct('orig',[], 'vec',[0,0,1]');
% for i=1:V.dim(1)
%     for j=1:V.dim(2)
%         ray.orig = [i,j,1]';
%         [~, P] = spm_mesh_ray_intersect(v, ray);
%         h(i,j,1:size(P,1)) = P(:,3);
%     end
% end

[x,y] = ndgrid(1:V.dim(1),1:V.dim(2));
ray = [x(:)';y(:)';repmat([0,0,0,1]',1,numel(x))];
[~,h] = spm_mesh_ray_triangle(double(v), ray);
h = reshape(h,V.dim(1),V.dim(2),[]);

%-Count number of intersections at Z > z
%--------------------------------------------------------------------------
R = false(V.dim);

for k=1:V.dim(3)
    R(:,:,k) = mod(sum(h > k, 3), 2);
end


%==========================================================================
function R = voxelise_with_inside_mesh_test(M, V)
I = 1:prod(V.dim);
[X,Y,Z] = ind2sub(V.dim, I');
XYZ = spm_mesh_transform([X,Y,Z], V.mat);

R = false(V.dim);

for i=1:size(XYZ,1)
    R(i) = spm_mesh_inside(M,XYZ(i,:));
    fprintf('%d out of %d\n',i,size(XYZ,1));
end
