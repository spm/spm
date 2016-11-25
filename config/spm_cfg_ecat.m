function ecat = spm_cfg_ecat
% SPM Configuration file for ECAT Import
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_ecat.m 6952 2016-11-25 16:03:13Z guillaume $


%--------------------------------------------------------------------------
% data ECAT files
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'ECAT files';
data.help    = {'Select the ECAT files to convert.'};
data.filter  = 'any';
data.ufilter = '.*v';
data.num     = [1 Inf];

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
ext.values = {'img', 'nii'};
ext.def    = @(val)spm_get_defaults('images.format', val{:});

%--------------------------------------------------------------------------
% opts Options
%--------------------------------------------------------------------------
opts         = cfg_branch;
opts.tag     = 'opts';
opts.name    = 'Options';
opts.val     = {ext };
opts.help    = {'Conversion options'};

%--------------------------------------------------------------------------
% ecat ECAT Import
%--------------------------------------------------------------------------
ecat         = cfg_exbranch;
ecat.tag     = 'ecat';
ecat.name    = 'ECAT Import';
ecat.val     = {data opts };
ecat.help    = {
    'ECAT 7 Conversion.'
    'ECAT 7 is the image data format used by the more recent CTI PET scanners.'
    }';
ecat.prog    = @convert_ecat;
ecat.vout    = @vout;
ecat.modality = {'PET'};


%==========================================================================
function out = convert_ecat(job)
for i=1:numel(job.data)
    spm_ecat2nifti(job.data{i},job.opts);
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
