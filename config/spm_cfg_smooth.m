function smooth = spm_cfg_smooth
% SPM Configuration file for Smooth
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_smooth.m 6148 2014-09-03 15:49:04Z guillaume $


%--------------------------------------------------------------------------
% data Images to Smooth
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Images to Smooth';
data.help    = {'Specify the images to smooth. The smoothed images are written to the same subdirectories as the original images and are prefixed with a ''s''. The prefix can be changed by an option setting.'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [0 Inf];
data.preview = @(f) spm_check_registration(char(f));

%--------------------------------------------------------------------------
% fwhm FWHM
%--------------------------------------------------------------------------
fwhm         = cfg_entry;
fwhm.tag     = 'fwhm';
fwhm.name    = 'FWHM';
fwhm.help    = {'Specify the full-width at half maximum (FWHM) of the Gaussian smoothing kernel in mm. Three values should be entered, denoting the FWHM in the x, y and z directions.'};
fwhm.strtype = 'r';
fwhm.num     = [1 3];
fwhm.def     = @(val)spm_get_defaults('smooth.fwhm', val{:});

%--------------------------------------------------------------------------
% dtype Data Type
%--------------------------------------------------------------------------
dtype         = cfg_menu;
dtype.tag     = 'dtype';
dtype.name    = 'Data Type';
dtype.help    = {'Data-type of output images.  SAME indicates the same datatype as the original images.'};
dtype.labels  = {
                'SAME'
                'UINT8   - unsigned char'
                'INT16   - signed short'
                'INT32   - signed int'
                'FLOAT32 - single prec. float'
                'FLOAT64 - double prec. float'
}';
dtype.values  = {0 spm_type('uint8') spm_type('int16') spm_type('int32') spm_type('float32') spm_type('float64')};
dtype.val     = {0};

%--------------------------------------------------------------------------
% im Implicit masking
%--------------------------------------------------------------------------
im         = cfg_menu;
im.tag     = 'im';
im.name    = 'Implicit masking';
im.help    = {'An "implicit mask" is a mask implied by a particular voxel value (0 for images with integer type, NaN for float images).'
              'If set to ''Yes'', the implicit masking of the input image is preserved in the smoothed image.'};
im.labels  = {'Yes' 'No'};
im.values  = {1 0};
im.val     = {0};

%--------------------------------------------------------------------------
% prefix Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the smoothed image file(s). Default prefix is ''s''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('smooth.prefix', val{:});

%--------------------------------------------------------------------------
% smooth Smooth
%--------------------------------------------------------------------------
smooth       = cfg_exbranch;
smooth.tag   = 'smooth';
smooth.name  = 'Smooth';
smooth.val   = {data fwhm dtype im prefix};
smooth.help  = {'This is for smoothing (or convolving) image volumes with a Gaussian kernel of a specified width. It is used as a preprocessing step to suppress noise and effects due to residual differences in functional and gyral anatomy during inter-subject averaging.'};
smooth.prog  = @spm_run_smooth;
smooth.vout  = @vout;


%==========================================================================
function dep = vout(varargin)
% Output file names will be saved in a struct with field .files
dep(1)            = cfg_dep;
dep(1).sname      = 'Smoothed Images';
dep(1).src_output = substruct('.','files');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
