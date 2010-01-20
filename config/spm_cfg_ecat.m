function ecat = spm_cfg_ecat
% SPM Configuration file for ECAT Import
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_ecat.m 3691 2010-01-20 17:08:30Z guillaume $

rev = '$Rev: 3691 $';
% ---------------------------------------------------------------------
% data ECAT files
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'ECAT files';
data.help    = {'Select the ECAT files to convert.'};
data.filter = 'any';
data.ufilter = '.*v';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% ext Output image format
% ---------------------------------------------------------------------
ext         = cfg_menu;
ext.tag     = 'ext';
ext.name    = 'Output image format';
ext.help    = {'Output files can be written as .img + .hdr, or the two can be combined into a .nii file.'};
ext.labels = {
              'Two file (img+hdr) NIfTI'
              'Single file (nii) NIfTI'
}';
ext.values = {
              'img'
              'nii'
}';
ext.def    = @(val)spm_get_defaults('images.format', val{:});
% ---------------------------------------------------------------------
% opts Options
% ---------------------------------------------------------------------
opts         = cfg_branch;
opts.tag     = 'opts';
opts.name    = 'Options';
opts.val     = {ext };
opts.help    = {'Conversion options'};
% ---------------------------------------------------------------------
% ecat ECAT Import
% ---------------------------------------------------------------------
ecat         = cfg_exbranch;
ecat.tag     = 'ecat';
ecat.name    = 'ECAT Import';
ecat.val     = {data opts };
ecat.help    = {'ECAT 7 Conversion.  ECAT 7 is the image data format used by the more recent CTI PET scanners.'};
ecat.prog = @convert_ecat;
ecat.modality = {'PET'};
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function convert_ecat(job)
for i=1:length(job.data),
    spm_ecat2nifti(job.data{i},job.opts);
end;
return;
