function opts = spm_config_delete
% Configuration file for deleting files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_config_delete.m 1150 2008-02-14 16:16:42Z guillaume $

deletefiles.type   = 'files';
deletefiles.name   = 'Files to delete';
deletefiles.tag    = 'deletefiles';
deletefiles.filter = '.*';
deletefiles.num    = [0 Inf];
deletefiles.help   = {'Select files to delete.'};

opts.type = 'branch';
opts.name = 'Delete Files';
opts.tag  = 'deletefiles';
opts.val  = {deletefiles};
opts.prog = @my_deletefiles;
opts.help = {[...
    'This facilty allows to delete files in a batch. Note that deleting ' ...
    'files will not make them disappear from file selection lists. ' ...
    'Therefore one has to be careful not to select the original files ' ...
    'after they have been programmed to be deleted.'],'',[...
    'If image files (.img or .nii) are selected, corresponding .hdr ' ...
    'or .mat files will be deleted as well, if they exist.']};

%------------------------------------------------------------------------
function my_deletefiles(varargin)
job = varargin{1};
for k = 1:numel(job.deletefiles)
    [p n e] = spm_fileparts(job.deletefiles{k});
    if strcmp(e,'.img') || strcmp(e,'.nii')
        spm_unlink(fullfile(p,[n '.hdr']));
        spm_unlink(fullfile(p,[n '.mat']));
    end
    spm_unlink(fullfile(p,[n e]));
end