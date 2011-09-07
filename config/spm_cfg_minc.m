function minc = spm_cfg_minc
% SPM Configuration file for MINC Import
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_minc.m 4466 2011-09-07 16:50:29Z guillaume $


%--------------------------------------------------------------------------
% data MINC files
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'MINC files';
data.help    = {'Select the MINC files to convert.'};
data.filter  = 'mnc';
data.ufilter = '.*';
data.num     = [1 Inf];

%--------------------------------------------------------------------------
% dtype Data Type
%--------------------------------------------------------------------------
dtype         = cfg_menu;
dtype.tag     = 'dtype';
dtype.name    = 'Data Type';
dtype.help    = {'Data-type of output images. Note that the number of bits used determines the accuracy, and the amount of disk space needed.'};
dtype.labels  = {
                'UINT8   - unsigned char'
                'INT16   - signed short'
                'INT32   - signed int'
                'FLOAT32 - single prec. float'
                'FLOAT64 - double prec. float'
}';
dtype.values = {spm_type('uint8') spm_type('int16') spm_type('int32') ...
                spm_type('float32') spm_type('float64')};
dtype.val    = {spm_type('int16')};

%--------------------------------------------------------------------------
% ext Output image format
%--------------------------------------------------------------------------
ext         = cfg_menu;
ext.tag     = 'ext';
ext.name    = 'Output image format';
ext.help    = {'Output files can be written as .img + .hdr, or the two can be combined into a .nii file.'};
ext.labels = {
              'Two file (img+hdr) NIfTI'
              'Single file (nii) NIfTI'
}';
ext.values = { 'img', 'nii' };
ext.def    = @(val)spm_get_defaults('images.format', val{:});

%--------------------------------------------------------------------------
% opts Options
%--------------------------------------------------------------------------
opts         = cfg_branch;
opts.tag     = 'opts';
opts.name    = 'Options';
opts.val     = {dtype ext};
opts.help    = {'Conversion options'};

%--------------------------------------------------------------------------
% minc MINC Import
%--------------------------------------------------------------------------
minc         = cfg_exbranch;
minc.tag     = 'minc';
minc.name    = 'MINC Import';
minc.val     = {data opts};
minc.help    = {'MINC Conversion.  MINC is the image data format used for exchanging data within the ICBM community, and the format used by the MNI software tools. It is based on NetCDF. MINC is no longer supported for reading images into SPM, so MINC files need to be converted to NIFTI format in order to use them. See http://www.bic.mni.mcgill.ca/software/ for more information.'};
minc.prog    = @spm_run_minc;
minc.vout    = @vout;


%==========================================================================
function out = spm_run_minc(job)
for i=1:numel(job.data)
    spm_mnc2nifti(job.data{i},job.opts);
end

out.files = cell(size(job.data));
for i=1:numel(job.data)
    out.files{i} = spm_file(job.data{i}, 'path',pwd, 'ext',job.opts.ext);
end


%==========================================================================
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'Converted Images';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
