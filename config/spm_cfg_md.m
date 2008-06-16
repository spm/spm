function md = spm_cfg_md
% SPM Configuration file for making directory function
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_md.m 1827 2008-06-16 13:54:37Z guillaume $

% ----------------------------------------------------------------------
% basedir Select a base directory
% ----------------------------------------------------------------------
basedir         = cfg_files;
basedir.tag     = 'basedir';
basedir.name    = 'Select a base directory';
basedir.help    = {'Select a base directory.'};
basedir.filter  = 'dir';
basedir.ufilter = '.*';
basedir.num     = [1 1];

% ----------------------------------------------------------------------
% name Enter a directory name
% ----------------------------------------------------------------------
name            = cfg_entry;
name.tag        = 'name';
name.name       = 'Enter a directory name';
name.help       = {'Enter a directory name'};
name.strtype    = 's';
name.num        = [1 Inf];

% ----------------------------------------------------------------------
% md Make Directory (Deprecated)
% ----------------------------------------------------------------------
md              = cfg_exbranch;
md.tag          = 'md';
md.name         = 'Make Directory (Deprecated)';
md.val          = {basedir name};
md.help         = {
    'This facilty allows programming a directory change. Directories are selected in the right listbox.'
    'This module is deprecated and has been moved to BasicIO.'
    'Jobs which are ready to run may continue using it, but the module inputs can not be changed via GUI. Please switch to the BasicIO module instead.'
}';
md.prog         = @my_mkdir;
md.hidden       = true;

%=======================================================================
function my_mkdir(varargin)
job = varargin{1};
if ~isempty(job.basedir) && ~isempty(job.name)
    mkdir(job.basedir{:},job.name);
end
