function cdir = spm_cfg_cdir
% SPM Configuration file for Change Directory
%__________________________________________________________________________

% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% -------------------------------------------------------------------------
% directory Select a directory
% -------------------------------------------------------------------------
directory         = cfg_files;
directory.tag     = 'directory';
directory.name    = 'Select a directory';
directory.help    = {'Select a directory to change to.'};
directory.filter  = 'dir';
directory.ufilter = '.*';
directory.num     = [1 1];

% -------------------------------------------------------------------------
% cdir Change Directory (Deprecated)
% -------------------------------------------------------------------------
cdir         = cfg_exbranch;
cdir.tag     = 'cdir';
cdir.name    = 'Change Directory (Deprecated)';
cdir.val     = {directory };
cdir.help    = {
                'This module is deprecated and has been moved to BasicIO.'
                'Jobs which are ready to run may continue using it, but the module inputs can not be changed via GUI. Please switch to the BasicIO module instead.'
                'This facilty allows programming a directory change. Directories are selected in the right listbox.'
}';
cdir.prog    = @my_job_cd;
cdir.hidden  = true;


%==========================================================================
function my_job_cd(varargin)
% job can be a job structure or the directory to change to.
job = varargin{1};
if isstruct(job)
    jobDir = job.directory;
else
    jobDir = job;
end
if ~isempty(jobDir)
    cd(char(jobDir));
    fprintf('New working directory: %s\n', char(jobDir));
end
