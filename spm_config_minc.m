function opts = spm_config_minc
% Configuration file for minc import jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_minc.m 1032 2007-12-20 14:45:55Z john $

%_______________________________________________________________________

data.type = 'files';
data.name = 'MINC files';
data.tag  = 'data';
data.filter = 'mnc';
data.num  = Inf;
data.help = {'Select the MINC files to convert.'};

dtype.type = 'menu';
dtype.name = 'Data Type';
dtype.tag  = 'dtype';
dtype.labels = {'UINT8  - unsigned char','INT16 - signed short','INT32 - signed int','FLOAT - single prec. float','DOUBLE - double prec. float'};
dtype.values = {spm_type('uint8'),spm_type('int16'),spm_type('int32'),spm_type('float32'),spm_type('float64')};
dtype.val = {spm_type('int16')};
dtype.help = {[...
'Data-type of output images. '...
'Note that the number of bits used determines '...
'the accuracy, and the amount of disk space needed.']};

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
op1.val  = {dtype,ext};
op1.help = {'Conversion options'};

opts.type = 'branch';
opts.name = 'MINC Import';
opts.tag  = 'minc';
opts.val  = {data,op1};
opts.prog = @convert_minc;
opts.vfiles = @vfiles;
opts.help = {[...
'MINC Conversion.  MINC is the image data format used for exchanging data '...
'within the ICBM community, and the format used by the MNI software tools. '...
'It is based on NetCDF, but due to be superceded by a new version relatively soon. '...
'MINC is no longer supported for reading images into SPM, so MINC files need to '...
'be converted to NIFTI format in order to use them. '...
'See http://www.bic.mni.mcgill.ca/software/ for more information.']};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function convert_minc(job)
for i=1:length(job.data),
    spm_mnc2nifti(job.data{i},job.opts);
end;
return;

function vf = vfiles(job)
vf  = cell(size(job.data));
for i=1:numel(job.data),
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    vf{i} = fullfile(pth,[nam job.opts.ext num]);
end;
