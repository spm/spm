function M = spm_mesh_transform(M,T,def)
% Apply a spatial transformation to vertices of a surface mesh
% FORMAT M = spm_mesh_transform(M,T,def)
% M        - a patch structure or a gifti object
% T        - a [4 x 4] transformation matrix [default: identity]
% def      - a deformation field (nifti object or filename) [default: none]
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_transform.m 4082 2010-10-07 16:11:48Z guillaume $

C = class(M);

if nargin < 2 || isempty(T)
    T = eye(4);
end

if nargin < 3
    M.vertices = [M.vertices ones(size(M.vertices,1),1)] * T(1:3,:)';
else
    M = spm_swarp(M,def,T);
end

M = feval(C, M);
