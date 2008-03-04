function out = spm_run_dicom(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_dicom.m 1185 2008-03-04 16:31:21Z volkmar $


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
out = spm_dicom_convert(hdr,'all',root_dir,job.convopts.format);

if ~isempty(job.outdir)
    fprintf('   Changing back to directory: %s\n', wd);
    cd(wd);
end

