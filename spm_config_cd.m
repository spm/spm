function opts = spm_config_cd
% Configuration file for changing directory function
%_______________________________________________________________________
% Darren Gitelmann $Id$

%_______________________________________________________________________

w = spm_jobman('HelpWidth');

%_______________________________________________________________________

data.type = 'dirs';
data.name = 'Select a directory';
data.tag  = 'directory';
data.filter = '.';
data.num  = 1;
data.help = {'Select a directory to change to.'};

opts.type = 'branch';
opts.name = 'Change Directory';
opts.tag  = 'cd';
opts.val  = {data};
opts.prog = @my_cd;
p1 = spm_justify(w,[...
'This facilty allows programming a directory change. Directories are ',...
'selected in the right listbox.']);

opts.help = {p1{:}};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_cd(varargin)
job = varargin{1};
if ~isempty(job.directory)
    cd(job.directory{:});
end
return;
