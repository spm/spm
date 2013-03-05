function bbox = spm_cfg_bbox
% SPM Configuration file for Get Bounding Box
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_cfg_bbox.m 5301 2013-03-05 18:33:39Z ged $

img             = cfg_files;
img.tag         = 'image';
img.name        = 'Image';
img.help        = {'Image for which to determine bounding box.'};
img.filter      = 'image';
img.ufilter     = '.*';
img.num         = [1 1];

fov             = cfg_const;
fov.tag         = 'fov';
fov.name        = 'Field of view';
fov.help        = {'Bounding box is for entire field of view of image.'};
fov.val         = {'fv'};

thr             = cfg_entry;
thr.tag         = 'threshold';
thr.name        = 'Threshold';
thr.help        = {'Bounding box is for set of voxels above this value.'};
thr.strtype     = 'r';
thr.num         = [1 1];

stv             = cfg_branch;
stv.tag         = 'stv';
stv.name        = 'Supra-threshold voxels';
stv.help        = {'Bounding box is for voxels over specified threshold.'};
stv.val         = {thr};

bbdef           = cfg_choice;
bbdef.tag       = 'bbdef';
bbdef.name      = 'Bounding box definition';
bbdef.help      = {
    ['Bounding box for field of view (using only header information) ' ...
    'or for the set of voxels over a specified threshold.']};
bbdef.values    = {fov stv};
bbdef.val       = {fov};

bbox            = cfg_exbranch;
bbox.tag        = 'bbox';
bbox.name       = 'Get Bounding Box';
bbox.val        = {img bbdef};
bbox.help       = {
    ['Determine the bounding box of an image, i.e. the [2 x 3] array ' ...
    'of the minimum and maximum X, Y, and Z coordinates (in mm), ']
    'BB = [min_X min_Y min_Z'
    '      max_X max_Y max_Z]'};
bbox.prog       = @(job) bbox_run(job);
bbox.vout       = @(job) bbox_vout(job);


%--------------------------------------------------------------------------
function out = bbox_run(job)
try
    thr = job.bbdef.stv.threshold;
catch
    thr = 'fv';
end

out.bb = spm_get_bbox(char(job.image), thr);

fprintf('\nImage:\n\t%s\n\nBounding box:\n', char(job.image));
disp(out.bb)

%--------------------------------------------------------------------------
function vout = bbox_vout(job)
vout            = cfg_dep;
vout.sname      = 'BB';
vout.src_output = substruct('.', 'bb');
vout.tgt_spec   = cfg_findspec({{'strtype','e'}});
