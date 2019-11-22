function voi = bf_sources_voi(BF, S)
% Generate a set of VOIs specified in MNI coordinates
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% $Id: bf_sources_voi.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the VOI'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates';
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the VOI in MNI coordinates'};
    pos.val = {};
    
    ori = cfg_entry;
    ori.tag = 'ori';
    ori.name = 'Orientation';
    ori.strtype = 'r';
    ori.num = [1 3];
    ori.help = {'Source orientatons (only for single points, leave zeros for unoriented)'};
    ori.val = {[0 0 0]};
    
    voidef = cfg_branch;
    voidef.tag = 'voidef';
    voidef.name = 'VOI';
    voidef.val = {label, pos, ori};
    
    mask = cfg_files;
    mask.tag = 'mask';
    mask.name = 'MNI mask';
    mask.filter = 'image';
    mask.ufilter = '.*';
    mask.num     = [1 1];
    mask.help = {'Select a mask image'};
    
    maskdef = cfg_branch;
    maskdef.tag = 'maskdef';
    maskdef.name = 'Mask VOI';
    maskdef.val  = {label, mask};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'VOIs';
    vois.num  = [1 Inf];
    vois.values = {voidef, maskdef};
    vois.val = {voidef};
    
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the VOIs (leave 0 for single point)'};
    
    resolution = cfg_entry;
    resolution.tag = 'resolution';
    resolution.name = 'Resolution';
    resolution.strtype = 'r';
    resolution.num = [1 1];
    resolution.val = {5};
    resolution.help = {'Resolution for placing grid points in each VOI (in mm)'};
    
    voi = cfg_branch;
    voi.tag = 'voi';
    voi.name = 'VOIs in MNI space';
    voi.val = {vois, radius, resolution};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end


iskull = export(gifti(BF.data.mesh.tess_iskull), 'ft');

M1 = BF.data.transforms.toNative;
M1 = BF.data.transforms.toMNI/M1;

iskull = ft_convert_units(ft_transform_geometry(M1, iskull));

       
% transform MNI coords in MNI space into space where we are doing the
% beamforming
M = inv(BF.data.transforms.toMNI);

if S.radius > 0
    vec = -S.radius:S.resolution:S.radius;
    [X, Y, Z]  = ndgrid(vec, vec, vec);
    sphere   = [X(:) Y(:) Z(:)];
    sphere(sqrt(X(:).^2 + Y(:).^2 + Z(:).^2) > S.radius, :) = [];
    npnt = size(sphere, 1);
else
    sphere = 0;
    npnt = 1;
end

grid = bf_sources_grid(BF, struct('resolution', S.resolution, 'space', 'MNI template'));
mnigrid = ft_transform_geometry(BF.data.transforms.toMNI, grid);

nvoi = numel(S.vois);
voi = [];
voi.label = {};
voi.pos = [];
ori = [];
voi.pos2voi = [];

for i = 1:nvoi
    switch char(fieldnames(S.vois{i}))
        case 'voidef'
            voi.label{i} = S.vois{i}.voidef.label;
            voi.pos = [voi.pos; sphere+repmat(S.vois{i}.voidef.pos, npnt, 1)];
            ori     = [ori; S.vois{i}.voidef.ori];
            voi.pos2voi  = [voi.pos2voi i*ones(1, npnt)];
        case 'maskdef'
            voi.label{i} = S.vois{i}.maskdef.label;
            V   = spm_vol(char(S.vois{i}.maskdef.mask));
            
            vox = spm_eeg_inv_transform_points(inv(V.mat), mnigrid.pos);
            Y   = spm_sample_vol(V, vox(:, 1),  vox(:, 2), vox(:, 3), 0);
            ind = find(~isnan(Y) & abs(Y)>0);
            voi.pos = [voi.pos; mnigrid.pos(ind, :)];
            ori = [ori;zeros(length(ind), 3)];
            voi.pos2voi  = [voi.pos2voi i*ones(1, length(ind))];
    end
end

voi.label = voi.label(:);

% Remove points outside the brain
inside = ft_inside_headmodel(voi.pos, struct('bnd', iskull));

voi.pos(~inside, :)  = [];
voi.pos2voi(~inside) = [];

if any(any(ori))
    voi.ori = ori(inside, :);
end


voi = ft_transform_geometry(M, voi);