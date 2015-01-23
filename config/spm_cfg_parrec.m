function parrec = spm_cfg_parrec
% SPM Configuration file for Philips PAR/REC Import
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_cfg_parrec.m 6315 2015-01-23 17:09:06Z guillaume $


%--------------------------------------------------------------------------
% data PAR files
%--------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'PAR files';
data.help    = {'Select the PAR files to convert.'};
data.filter  = 'any';
data.ufilter = '.*(par|PAR)';
data.num     = [1 Inf];

%--------------------------------------------------------------------------
% outdir Output directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where files are written.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

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
opts.val     = {outdir ext};
opts.help    = {'Conversion options'};

%--------------------------------------------------------------------------
% parrec ECAT Import
%--------------------------------------------------------------------------
parrec         = cfg_exbranch;
parrec.tag     = 'parrec';
parrec.name    = 'PAR/REC Import';
parrec.val     = {data opts};
parrec.help    = {'Philips PAR/REC Import.'};
parrec.prog    = @convert_parrec;
parrec.vout    = @vout;


%==========================================================================
function out = convert_parrec(job)
for i=1:numel(job.data)
    spm_parrec2nifti(job.data{i},job.opts);
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
