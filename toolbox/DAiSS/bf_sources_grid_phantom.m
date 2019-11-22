function res = bf_sources_grid_phantom(BF, S)
% Generate beamforming grid
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_sources_grid_phantom.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0 
    resolution = cfg_entry;
    resolution.tag = 'resolution';
    resolution.name = 'Grid resolution';
    resolution.strtype = 'n';
    resolution.num = [1 1];
    resolution.val = {5};
    resolution.help = {'Select the resolution of the grid (in mm)'};       

    grid = cfg_branch;
    grid.tag = 'grid_phantom';
    grid.name = 'Phantom Grid';
    grid.val = {resolution};
    
    res = grid;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

vol = BF.data.MEG.vol;

if ~isequal(vol.type,  'singlesphere') || any(vol.o)
    error('Expecting a spherical volume centered at the origin');
end

M1 = BF.data.transforms.toNative;
M1 = BF.data.transforms.toHead/M1;
M2 = inv(BF.data.transforms.toHead);

mn = -vol.r*[1 1 1];
mx =  vol.r*[1 1 1];

resolution = S.resolution;

if isequal(vol.unit, 'm')
    resolution = 1e-3*resolution;
end

% If zero is inside the brain, make sure grid points fall on multiples of
% resolution to ease simulating data from points on the grid
if mn(1)<0 && mx(1)>0
    grid.xgrid = [fliplr(0:-resolution:mn(1)) resolution:resolution:mx(1)];
else
    grid.xgrid = mn(1):resolution:mx(1);
end

if mn(2)<0 && mx(2)>0
    grid.ygrid = [fliplr(0:-resolution:mn(2)) resolution:resolution:mx(2)];
else
    grid.ygrid = mn(2):resolution:mx(2);
end

if mn(3)<0 && mx(3)>0
    grid.zgrid = [fliplr(0:-resolution:mn(3)) resolution:resolution:mx(3)];
else
    grid.zgrid = mn(3):resolution:mx(3);
end

grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
[X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);

pos   = [X(:) Y(:) Z(:)];

inside = ft_inside_headmodel(pos, vol);

pos    = spm_eeg_inv_transform_points(M2, pos);

grid.allpos  = pos;
grid.inside  = find(inside);
grid.outside = find(~inside);

grid.pos     = pos(inside, :);

res = grid;