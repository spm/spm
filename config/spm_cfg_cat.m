function cat = spm_cfg_cat
% SPM Configuration file for 3D to 4D volumes conversion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_cfg_cat.m 5828 2014-01-03 18:38:35Z guillaume $

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
cat      = cfg_exbranch;
cat.tag  = 'cat';
cat.name = '3D to 4D File Conversion';
cat.val  = {vols name dtype};
cat.help = {'Concatenate a number of 3D volumes into a single 4D file.'};
cat.prog = @(job)spm_run_cat('run',job);
cat.vout = @(job)spm_run_cat('vout',job);


%==========================================================================
function out = spm_run_cat(cmd, job)

switch lower(cmd)
    case 'run'
        V                 = char(job.vols{:});
        dt                = job.dtype;
        fname             = job.name;
        V4                = spm_file_merge(V,fname,dt);
        out.mergedfile    = {V4(1).fname};
        
    case 'vout'
        out(1)            = cfg_dep;
        out(1).sname      = 'Concatenated 4D Volume';
        out(1).src_output = substruct('.','mergedfile');
        out(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
