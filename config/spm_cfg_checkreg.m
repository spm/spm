function checkreg = spm_cfg_checkreg
% SPM Configuration file for Check Reg
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_checkreg.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% data Images to Display
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Images to Display';
data.help    = {'Images to display.'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 24];

%--------------------------------------------------------------------------
% checkreg Check Registration
%--------------------------------------------------------------------------
checkreg      = cfg_exbranch;
checkreg.tag  = 'checkreg';
checkreg.name = 'Check Registration';
checkreg.val  = {data};
checkreg.help = {
    'Display of orthogonal views of one or more images.'
    ''
    'Clicking in any image moves the centre of the orthogonal views. Images are shown in orientations relative to that of the first selected image. The first specified image is shown at the top-left, and the last at the bottom right.  The fastest increment is in the left-to-right direction (the same as you are reading this).'
    ''
    'If you have put your images in the correct file format, then (possibly after specifying some rigid-body rotations):'
    '    The top-left image is coronal with the top (superior) of the head displayed at the top and the left shown on the left. This is as if the subject is viewed from behind.'
    '    The bottom-left image is axial with the front (anterior) of the head at the top and the left shown on the left. This is as if the subject is viewed from above.'
    '    The top-right image is sagittal with the front (anterior) of the head at the left and the top of the head shown at the top. This is as if the subject is viewed from the left.'
    }';
checkreg.prog = @check_reg;


%==========================================================================
function check_reg(job)
spm_check_registration(char(job.data));