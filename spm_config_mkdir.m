function opts = spm_config_mkdir
% Configuration file for making directory function
%_______________________________________________________________________
% DRG CS-RCS: $Id: spm_config_mkdir.m,v 1.1 2005-02-09 21:42:37-06 drg Exp drg $

%_______________________________________________________________________

w = spm_jobman('HelpWidth');

%_______________________________________________________________________

data.type = 'dirs';
data.name = 'Select a base directory';
data.tag  = 'basedir';
data.filter = '.';
data.num  = 1;
data.help = {'Select a directory to change to.'};

name.type  = 'entry';
name.name = 'Enter a directory name';
name.tag = 'name';
name.strtype = 's';
name.num     = [1 1];
name.help = 'Enter a directory name';

opts.type = 'branch';
opts.name = 'Make Directory';
opts.tag  = 'md';
opts.val  = {data,name};
opts.prog = @my_cd;
opts.vdirs = @vdirs_mydirs;
p1 = spm_justify(w,[...
'This facilty allows programming a directory change. Directories are ',...
'selected in the right listbox.']);

opts.help = {p1{:}};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_cd(varargin)
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

