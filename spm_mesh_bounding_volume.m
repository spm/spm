function bv = spm_mesh_bounding_volume(M,t)
% Bounding volume of a triangle mesh
% FORMAT bv = spm_mesh_bounding_volume(M,t)
% M        - a patch structure or GIfTI object
% t        - type of bounding volume [default: 'AABB']
%
% bv       - bounding volume
%__________________________________________________________________________
%
% See: https://en.wikipedia.org/wiki/Bounding_volume
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if nargin < 2
    t = 'AABB';
end

switch lower(t)
    case 'aabb'
        %-Axis-Aligned Bounding Box
        bv = double([min(M.vertices,[],1); max(M.vertices,[],1)]);
    case 'obb'
        %-Oriented Bounding Box
        error('Not implemented yet');
    case 'ellipsoid'
        %-Bounding ellipsoid
        error('Not implemented yet');
    case 'sphere'
        %-Bounding sphere
        error('Not implemented yet');
    case 'convhull'
        %-Convex hull
        bv = convhull(M.vertices(:,1),M.vertices(:,2),M.vertices(:,3));
    otherwise
        error('Unknown type of bounding volume.');
end
