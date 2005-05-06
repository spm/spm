function opts = spm_config_cd
% Configuration file for changing directory function
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Darren Gitelmann
% $Id: spm_config_cd.m 123 2005-05-06 12:15:13Z john $

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
opts.prog = @my_cd;
p1 = [...
'This facilty allows programming a directory change. Directories are ',...
'selected in the right listbox.'];

opts.help = {p1};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_cd(varargin)
job = varargin{1};
if ~isempty(job.directory),
    cd(job.directory{:});
end;
return;
