function opts = spm_config_cd
% Configuration file for changing directory function
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman
% $Id: spm_config_cd.m 587 2006-08-07 04:38:22Z Darren $

data.type = 'files';
data.name = 'Select a directory';
data.tag  = 'directory';
data.filter = 'dir';
data.num  = 1;
data.help = {'Select a directory to change to.'};

opts.type = 'branch';
opts.name = 'Change Directory';
opts.tag  = 'cdir';
opts.val  = {data};
opts.prog = @my_job_cd;
p1 = [...
'This facilty allows programming a directory change. Directories are ',...
'selected in the right listbox.'];

opts.help = {p1};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_job_cd(varargin)
% job can be a job structure or the directory to change to.
job = varargin{1};
if isstruct(job)
    jobDir = job.directory;
else
    jobDir = job;
end
if ~isempty(jobDir),
    cd(char(jobDir));
    fprintf('New working directory: %s\n', char(jobDir));
end;
return;
