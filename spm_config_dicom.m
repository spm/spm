function opts = spm_config_dicom
% Configuration file for dicom import jobs
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_dicom.m 506 2006-04-27 14:46:29Z volkmar $


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
'Select a directory where files are written.'...
'Default is current directory.']};

opts.type = 'branch';
opts.name = 'DICOM Import';
opts.tag  = 'dicom';
opts.val  = {data,outdir};
opts.prog = @convert_dicom;
opts.help = {[...
'DICOM Conversion.  Most scanners produce data in DICOM format. '...
'This routine attempts to convert DICOM files into SPM compatible '...
'image volumes, which are written into the current directory by '...
'default. Note that not all flavours of DICOM can be handled, as '...
'DICOM is a very complicated format, and some scanner manufacturers '...
'use their own fields, which are not in the official documentation '...
'at http://medical.nema.org/']};

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

