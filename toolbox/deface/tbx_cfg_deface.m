function job = tbx_cfg_deface
% SPM Configuration file for toolbox 'De-Face'
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: tbx_cfg_deface.m 5713 2013-10-29 15:23:18Z john $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','deface')); end

images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images to de-face';
images.help    = {'Specify the NIfTI images to strip the face from.'};
images.filter = 'nifti';
images.ufilter = '.*';
images.num     = [0 Inf];

% Select the images
job         = cfg_exbranch;
job.tag     = 'deface';
job.name    = 'De-face Images';
job.val     = {images};
job.help    = {'This is a little routine for attempting to strip the face from images, so individuals are more difficult to identify from surface renderings.'};
job.prog = @spm_deface;
job.vout = @vout;

function dep = vout(varargin)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'De-faced images';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

