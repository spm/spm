function M = spm_mesh_transform(M,T,def)
% Apply a spatial transformation to vertices of a surface mesh
% FORMAT M = spm_mesh_transform(M,T,def)
% M        - a patch structure or a gifti object or [nv x 3] array
% T        - a [4 x 4] transformation matrix [default: identity]
% def      - a deformation field (nifti object or filename) [default: none]
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


C = class(M);

if nargin < 2 || isempty(T)
    T = eye(4);
end

if nargin < 3
    if isfield(M,'vertices')
        M.vertices = [M.vertices ones(size(M.vertices,1),1)] * T(1:3,:)';
    else
        M = [M ones(size(M,1),1)] * T(1:3,:)';
    end
else
    M = spm_swarp(M,def,T);
end

M = feval(C, M);
