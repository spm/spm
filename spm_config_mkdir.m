function opts = spm_config_mkdir
% Configuration file for making directory function
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelman
% $Id: spm_config_mkdir.m 472 2006-03-08 17:48:52Z john $

%_______________________________________________________________________


data.type = 'files';
data.name = 'Select a base directory';
data.tag  = 'basedir';
data.filter = 'dir';
data.num  = 1;
data.help = {'Select a base directory.'};

name.type  = 'entry';
name.name = 'Enter a directory name';
name.tag = 'name';
name.strtype = 's';
name.num     = [1 1];
name.help = {'Enter a directory name'};

opts.type = 'branch';
opts.name = 'Make Directory';
opts.tag  = 'md';
opts.val  = {data,name};
opts.prog = @my_mkdir;
opts.vfiles = @vdirs_mydirs;
opts.help = {[...
'This facilty allows programming a directory change. Directories are ',...
'selected in the right listbox.']};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_mkdir(varargin)
job = varargin{1};
if ~isempty(job.basedir) && ~isempty(job.name)
    mkdir(job.basedir{:},job.name);
end
%------------------------------------------------------------------------

%------------------------------------------------------------------------

function vd = vdirs_mydirs(varargin)

job = varargin{1};

vd = {fullfile(job.basedir{:},[job.name,filesep])};
return;

