function opts = spm_config_ecat
% Configuration file for ecat import jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_ecat.m 512 2006-05-05 08:14:50Z volkmar $

%_______________________________________________________________________

data.type = 'files';
data.name = 'ECAT files';
data.tag  = 'data';
data.ufilter = '.*v';
data.filter = 'any';
data.num  = Inf;
data.help = {'Select the ECAT files to convert.'};

dtype.type = 'menu';
dtype.name = 'Data Type';
dtype.tag  = 'dtype';
dtype.labels = {'UINT8  - unsigned char','INT16 - signed short','INT32 - signed int','FLOAT - single prec. float','DOUBLE - double prec. float'};
dtype.values = {spm_type('uint8'),spm_type('int16'),spm_type('int32'),spm_type('float32'),spm_type('float64')};
dtype.val = {spm_type('int16')};
dtype.help = {[...
'Data-type of output images. '...
'Note that the number of bits used determines '...
'the accuracy, and the abount of disk space needed.']};

ext.type = 'menu';
ext.name = 'NIFTI Type';
ext.tag  = 'ext';
ext.labels = {'.nii only','.img + .hdr'};
ext.values = {'.nii','.img'};
ext.val    = {'.img'};
ext.help   = {[...
'Output files can be written as .img + .hdr, ',...
'or the two can be combined into a .nii file.']};

op1.type = 'branch';
op1.name = 'Options';
op1.tag  = 'opts';
op1.val  = {ext};
op1.help = {'Conversion options'};

opts.type = 'branch';
opts.name = 'ECAT Import';
opts.tag  = 'ecat';
opts.val  = {data,op1};
opts.prog = @convert_ecat;
opts.modality = {'PET'};
%opts.vfiles = @vfiles;
opts.help = {[...
'ECAT 7 Conversion.  ECAT 7 is the image data format used by the more recent CTI '...
'PET scanners.']};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function convert_ecat(job)
for i=1:length(job.data),
    spm_ecat2nifti(job.data{i},job.opts);
end;
return;

%function vf = vfiles(job)
%vf  = cell(size(job.data));
%for i=1:numel(job.data),
%    [pth,nam,ext,versn] = fileparts(job.data{i});
%    vf{i} = fullfile(pwd,[nam job.opts.ext versn]);
%end;
