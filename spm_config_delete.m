function opts = spm_config_delete
% Configuration file for deleting files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_config_delete.m 1115 2008-01-22 15:58:01Z guillaume $

deletefiles.type = 'files';
deletefiles.name = 'Files to delete';
deletefiles.tag  = 'deletefiles';
deletefiles.filter = '.*';
deletefiles.num  = [0 Inf];
deletefiles.help = {'Select files to delete.'};

opts.type = 'branch';
opts.name = 'Delete Files';
opts.tag  = 'deletefiles';
opts.val  = {deletefiles};
opts.prog = @my_deletefiles;
opts.vfiles = @vfiles_deletefiles;
opts.help = {[...
    'This facilty allows to delete files in a batch.']};

%------------------------------------------------------------------------
function my_deletefiles(varargin)
job = varargin{1};
for k = 1:numel(job.deletefiles)
    spm_unlink(job.deletefiles{k});
end