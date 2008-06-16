function movefile = spm_cfg_movefile
% SPM Configuration file for move file function
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_movefile.m 1827 2008-06-16 13:54:37Z guillaume $

% ----------------------------------------------------------------------
% srcfiles Files to move
% ----------------------------------------------------------------------
srcfiles         = cfg_files;
srcfiles.tag     = 'srcfiles';
srcfiles.name    = 'Files to move';
srcfiles.help    = {'Select files to move.'};
srcfiles.filter  = '.*';
srcfiles.ufilter = '.*';
srcfiles.num     = [0 Inf];

% ----------------------------------------------------------------------
% targetdir Target directory
% ----------------------------------------------------------------------
targetdir         = cfg_files;
targetdir.tag     = 'targetdir';
targetdir.name    = 'Target directory';
targetdir.help    = {'Select target directory.'};
targetdir.filter  = 'dir';
targetdir.ufilter = '.*';
targetdir.num     = [1 1];

% ----------------------------------------------------------------------
% movefile Move Files (Deprecated)
% ----------------------------------------------------------------------
movefile          = cfg_exbranch;
movefile.tag      = 'movefile';
movefile.name     = 'Move Files (Deprecated)';
movefile.val      = {srcfiles targetdir };
movefile.help     = {
    'This module is deprecated and has been moved to BasicIO.'
    'Jobs which are ready to run may continue using it, but the module inputs can not be changed via GUI. Please switch to the BasicIO module instead.'
    ''
    'This facilty allows to move files in a batch. Note that moving files will not make them disappear from file selection lists. Therefore one has to be careful not to select the original files after they have been programmed to be moved.'
    ''
    'If image files (.*img or .*nii) are selected, corresponding hdr or mat files will be moved as well, if they exist.'
};
movefile.prog     = @my_movefile;
movefile.hidden   = true;

%=======================================================================
function my_movefile(varargin)
job = varargin{1};
for k = 1:numel(job.srcfiles)
    [p n e v] = spm_fileparts(job.srcfiles{k});
    if strncmp(e,'.img',4)||strncmp(e,'.nii',4)
        try_movefile(fullfile(p,[n e]),job.targetdir{1});
        try_movefile(fullfile(p,[n '.mat']),job.targetdir{1});
        try_movefile(fullfile(p,[n '.hdr']),job.targetdir{1});
    else
        try_movefile(job.srcfiles{k},job.targetdir{1});
    end
end

%=======================================================================
function try_movefile(src,dest)
% silently try to move files
try
    movefile(src,dest);
end
