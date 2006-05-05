function opts = spm_config_dicom
% Configuration file for dicom import jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_dicom.m 511 2006-05-05 07:59:40Z volkmar $


%_______________________________________________________________________

data.type   = 'files';
data.name   = 'DICOM files';
data.tag    = 'data';
data.filter = '.*';
data.num    = Inf;
data.help   = {'Select the DICOM files to convert.'};

outdir.type   = 'files';
outdir.name   = 'Output directory';
outdir.tag    = 'outdir';
outdir.filter = 'dir';
outdir.num    = 1;
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

opts.type = 'branch';
opts.name = 'DICOM Import';
opts.tag  = 'dicom';
opts.val  = {data,root,outdir};
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

hdr = spm_dicom_headers(strvcat(job.data));
spm_dicom_convert(hdr);

if ~isempty(job.outdir)
    fprintf('   Changing back to directory: %s\n', wd);
    cd(wd);
end

