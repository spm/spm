function job = spm_cfg_deface
% SPM Configuration file for toolbox 'De-Face'
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_deface.m 6086 2014-07-03 16:08:44Z guillaume $


%--------------------------------------------------------------------------
% images Images to de-face
%--------------------------------------------------------------------------
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images to de-face';
images.help    = {'Specify the NIfTI images to strip the face from.'};
images.filter  = 'nifti';
images.ufilter = '.*';
images.num     = [0 Inf];

%--------------------------------------------------------------------------
% job De-face Images
%--------------------------------------------------------------------------
job       = cfg_exbranch;
job.tag   = 'deface';
job.name  = 'De-face Images';
job.val   = {images};
job.help  = {'This is a little routine for attempting to strip the face from images, so individuals are more difficult to identify from surface renderings.',...
             'De-faced images are prefixed by ''anon_''.'};
job.prog  = @spm_deface;
job.vout  = @vout;


%==========================================================================
function dep = vout(varargin)
% Output file names will be saved in a cell array of strings
dep(1)            = cfg_dep;
dep(1).sname      = 'De-faced images';
dep(1).src_output = substruct('()',{':'});
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
