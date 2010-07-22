function dicom = spm_cfg_dicom
% SPM Configuration file for DICOM Import
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_dicom.m 4012 2010-07-22 12:45:53Z volkmar $

% ---------------------------------------------------------------------
% data DICOM files
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'DICOM files';
data.help    = {'Select the DICOM files to convert.'};
data.filter  = 'any';
data.ufilter = '.*';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% root Directory structure for converted files
% ---------------------------------------------------------------------
root         = cfg_menu;
root.tag     = 'root';
root.name    = 'Directory structure for converted files';
root.help    = {'Choose root directory of converted file tree. The options are:'
                ''
                '* Output directory: ./<StudyDate-StudyTime>: Automatically determine the project name and try to convert into the output directory, starting with a StudyDate-StudyTime subdirectory. This option is useful if automatic project recognition fails and one wants to convert data into a project directory.'
                ''
                '* Output directory: ./<PatientID>: Convert into the output directory, starting with a PatientID subdirectory.'
                ''
                '* Output directory: ./<PatientName>: Convert into the output directory, starting with a PatientName subdirectory.'
                '* No directory hierarchy: Convert all files into the output directory, without sequence/series subdirectories'}';
root.labels  = {'Output directory: ./<StudyDate-StudyTime>/<ProtocolName>'
                'Output directory: ./<PatientID>/<ProtocolName>'
                'Output directory: ./<PatientID>/<StudyDate-StudyTime>/<ProtocolName>'
                'Output directory: ./<PatientName>/<ProtocolName>'
                'Output directory: ./<ProtocolName>'
                'No directory hierarchy'}';
root.values  = {'date_time'
                'patid'
                'patid_date'
                'patname'
                'series'
                'flat'}';
root.def     = @(val)spm_get_defaults('dicom.root', val{:});
% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where files are written.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];
% ---------------------------------------------------------------------
% format Output image format
% ---------------------------------------------------------------------
format         = cfg_menu;
format.tag     = 'format';
format.name    = 'Output image format';
format.help    = {'DICOM conversion can create separate img and hdr files or combine them in one file. The single file option will help you save space on your hard disk, but may be incompatible with programs that are not NIfTI-aware.'
                  'In any case, only 3D image files will be produced.'}';
format.labels  = {'Two file (img+hdr) NIfTI'
                  'Single file (nii) NIfTI'}';
format.values  = {'img' 'nii'};
format.def     = @(val)spm_get_defaults('images.format', val{:});
% ---------------------------------------------------------------------
% icedims Use ICEDims in filename
% ---------------------------------------------------------------------
icedims        = cfg_menu;
icedims.tag    = 'icedims';
icedims.name   = 'Use ICEDims in filename';
icedims.help   = {'If image sorting fails, one can try using the additional SIEMENS ICEDims information to create unique filenames. Use this only if there would be multiple volumes with exactly the same file names.'};
icedims.labels = {'No' 'Yes'};
icedims.values = {0 1};
icedims.val    = {0};
% ---------------------------------------------------------------------
% convopts Conversion options
% ---------------------------------------------------------------------
convopts       = cfg_branch;
convopts.tag   = 'convopts';
convopts.name  = 'Conversion options';
convopts.val   = {format icedims};
convopts.help  = {''};
% ---------------------------------------------------------------------
% dicom DICOM Import
% ---------------------------------------------------------------------
dicom          = cfg_exbranch;
dicom.tag      = 'dicom';
dicom.name     = 'DICOM Import';
dicom.val      = {data root outdir convopts};
dicom.help     = {'DICOM Conversion.  Most scanners produce data in DICOM format. This routine attempts to convert DICOM files into SPM compatible image volumes, which are written into the current directory by default. Note that not all flavours of DICOM can be handled, as DICOM is a very complicated format, and some scanner manufacturers use their own fields, which are not in the official documentation at http://medical.nema.org/'};
dicom.prog     = @spm_run_dicom;
dicom.vout     = @vout;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function dep = vout(job)
dep            = cfg_dep;
dep.sname      = 'Converted Images';
dep.src_output = substruct('.','files');
dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
