function res = bf_sources_grid(BF, S)
% Generate beamforming grid
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources_grid.m 4847 2012-08-16 17:29:23Z vladimir $

%--------------------------------------------------------------------------
if nargin == 0 
    resolution = cfg_entry;
    resolution.tag = 'resolution';
    resolution.name = 'Grid resolution';
    resolution.strtype = 'n';
    resolution.num = [1 1];
    resolution.val = {5};
    resolution.help = {'Select the resolution of the grid (in mm)'};
    
    space = cfg_menu;
    space.tag = 'space';
    space.name = 'Coordinate system';
    space.help = {'Select the coordinate system in which the grid should be generated'};
    space.labels = {'MNI template', 'MNI-aligned', 'Head', 'Native'};
    space.values = {'MNI template', 'MNI-aligned', 'Head', 'Native'};
    space.val = {'MNI template'};

    grid = cfg_branch;
    grid.tag = 'grid';
    grid.name = 'Grid';
    grid.val = {resolution, space};
    
    res = grid;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

iskull = export(gifti(BF.data.mesh.tess_iskull), 'ft');

M1 = BF.data.transforms.toMNI;

switch S.space
    case 'MNI template'
        M2 = inv(M1);
        M1 = eye(4);
    case 'MNI-aligned'
        M1 = BF.data.transforms.toMNI_aligned/M1;
        M2 = inv(BF.data.transforms.toMNI_aligned);
    case 'Head'
        M1 = BF.data.transforms.toHead/M1;
        M2 = inv(BF.data.transforms.toHead);
    case 'Native'
        M1 = BF.data.transforms.toNative/M1;
        M2 = inv(BF.data.transforms.toNative);
end        

iskull = ft_transform_geometry(M1, iskull);
mn = min(iskull.pnt);
mx = max(iskull.pnt);


grid.xgrid = mn(1):S.resolution:mx(1);
grid.ygrid = mn(2):S.resolution:mx(2);
grid.zgrid = mn(3):S.resolution:mx(3);

grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
[X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);

grid.pos   = [X(:) Y(:) Z(:)];

res = ft_transform_geometry(M2, grid);