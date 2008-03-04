function opts = spm_config_dicom
% Configuration file for dicom import jobs
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_config_dicom.m 1185 2008-03-04 16:31:21Z volkmar $


%_______________________________________________________________________

data.type   = 'files';
data.name   = 'DICOM files';
data.tag    = 'data';
data.filter = '.*';
data.num    = [1 Inf];
data.help   = {'Select the DICOM files to convert.'};

outdir.type   = 'files';
outdir.name   = 'Output directory';
outdir.tag    = 'outdir';
outdir.filter = 'dir';
outdir.num    = [1 1];
outdir.val    = {''};
outdir.help   = {[...
'Select a directory where files are written. '...
'Default is current directory.']};

root.type = 'menu';
root.name = 'Directory structure for converted files';
root.tag  = 'root';
root.labels = {'Output directory: ./<StudyDate-StudyTime>', ...
               'Output directory: ./<PatientID>', ...
               'Output directory: ./<PatientID>/<StudyDate-StudyTime>', ...
               'Output directory: ./<PatientName>', ...
               'No directory hierarchy'};
root.values = {'date_time',...
               'patid', 'patid_date', 'patname','flat'};
root.def  = 'util.dicom.root';
root.help = {['Choose root directory of converted file tree. The options ' ...
              'are:'], '',...
             ['* Output directory: ./<StudyDate-StudyTime>: ' ...
              'Automatically determine the project name and try to ' ...
              'convert into the output directory, starting with '... 
              'a StudyDate-StudyTime subdirectory. This option is useful if automatic '... 
              'project recognition fails and one wants to convert data into '...
              'a project directory.'], '',...
             ['* Output directory: ./<PatientID>: ' ...
              'Convert into the output directory, starting with '... 
              'a PatientID subdirectory.'], '',...
             ['* Output directory: ./<PatientName>: ' ...
              'Convert into the output directory, starting with '... 
              'a PatientName subdirectory.'],...
             ['* No directory hierarchy: Convert all files into the output ' ...
              'directory, without sequence/series subdirectories']};

format.type = 'menu';
format.name = 'Output image format';
format.tag  = 'format';
format.labels = {'Two file (img+hdr) NIfTI', 'Single file (nii) NIfTI'};
format.values = {'img', 'nii'};
format.val  = {'nii'};
format.help = {['DICOM conversion can create separate img and hdr files ' ...
                'or combine them in one file. The single file option will ' ...
                'help you save space on your hard disk, but may be ' ...
                'incompatible with programs that are not NIfTI-aware.'],...
               ['In any case, only 3D image files will be produced.']};

icedims.type = 'menu';
icedims.name = 'Use ICEDims in filename';
icedims.tag  = 'icedims';
icedims.labels = {'No','Yes'};
icedims.values = {0, 1};
icedims.val  = {0};
icedims.help = {['If image sorting fails, one can try using the additional ' ...
                'SIEMENS ICEDims information to create unique filenames. ' ...
                'Use this only if there would be multiple volumes with '...
                'exactly the same file names.']};

convopts.type = 'branch';
convopts.name = 'Conversion options';
convopts.tag  = 'convopts';
convopts.val  = {format,icedims};

opts.type = 'branch';
opts.name = 'DICOM Import';
opts.tag  = 'dicom';
opts.val  = {data,root,outdir,convopts};
opts.prog = @convert_dicom;
opts.help = {[...
'DICOM Conversion.  Most scanners produce data in DICOM format. '...
'This routine attempts to convert DICOM files into SPM compatible '...
'image volumes, which are written into the current directory by '...
'default. Note that not all flavours of DICOM can be handled, as '...
'DICOM is a very complicated format, and some scanner manufacturers '...
'use their own fields, which are not in the official documentation '...
'at http://medical.nema.org/']};

% add defaults, if not already set
global defaults
if isfield(defaults,'util')
  if isfield(defaults.util,'dicom')
    return;
  end;
end;
defaults.util.dicom.root = 'flat';

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function convert_dicom(job)

wd = pwd;
try
    if ~isempty(job.outdir{:})
        cd(job.outdir{:});
        fprintf('   Changing directory to: %s\n', job.outdir{:});
    end
catch
    error('Failed to change directory. Aborting DICOM import.');
end

if job.convopts.icedims
    root_dir = ['ice' job.root];
else
    root_dir = job.root;
end;

hdr = spm_dicom_headers(strvcat(job.data));
spm_dicom_convert(hdr,'all',root_dir,job.convopts.format);

if ~isempty(job.outdir)
    fprintf('   Changing back to directory: %s\n', wd);
    cd(wd);
end

