function cat = spm_cfg_cat
% SPM Configuration file for 3D to 4D volumes conversion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_cat.m 3613 2009-12-04 18:47:59Z guillaume $

%--------------------------------------------------------------------------
% vols 3D Volumes
%--------------------------------------------------------------------------
vols         = cfg_files;
vols.tag     = 'vols';
vols.name    = '3D Volumes';
vols.help    = {'Select the volumes to concatenate'};
vols.filter  = 'image';
vols.ufilter = '.*';
vols.num     = [1 Inf];

%--------------------------------------------------------------------------
% dtype Data Type
%--------------------------------------------------------------------------
dtype        = cfg_menu;
dtype.tag    = 'dtype';
dtype.name   = 'Data Type';
dtype.help   = {'Data-type of output image. SAME indicates the same datatype as the original images.'};
dtype.labels = {'SAME'
                'UINT8   - unsigned char'
                'INT16   - signed short'
                'INT32   - signed int'
                'FLOAT32 - single prec. float'
                'FLOAT64 - double prec. float'}';
dtype.values = {0 spm_type('uint8') spm_type('int16') spm_type('int32') spm_type('float32') spm_type('float64')};
dtype.val    = {spm_type('int16')}; % to match previous behaviour

%--------------------------------------------------------------------------
% name Output Filename
%--------------------------------------------------------------------------
name         = cfg_entry;
name.tag     = 'name';
name.name    = 'Output Filename';
name.help    = {'Specify the name of the output 4D volume file.'
                'A ''.nii'' extension will be added if not specified.'}';
name.strtype = 's';
name.num     = [1 Inf];
name.val     = {'4D.nii'};

%--------------------------------------------------------------------------
% cat 3D to 4D File Conversion
%--------------------------------------------------------------------------
cat         = cfg_exbranch;
cat.tag     = 'cat';
cat.name    = '3D to 4D File Conversion';
cat.val     = {vols name dtype};
cat.help    = {'Concatenate a number of 3D volumes into a single 4D file.'};
cat.prog = @spm_run_cat;
cat.vout = @vout;

%==========================================================================
function dep = vout(varargin)
% 4D output file will be saved in a struct with field .mergedfile
dep(1)            = cfg_dep;
dep(1).sname      = 'Concatenated 4D Volume';
dep(1).src_output = substruct('.','mergedfile');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
