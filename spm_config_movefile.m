function opts = spm_config_movefile
% Configuration file for move file function
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Volkmar Glauche
% $Id: spm_config_movefile.m 549 2006-06-07 12:37:29Z volkmar $

%_______________________________________________________________________


srcfiles.type = 'files';
srcfiles.name = 'Files to move';
srcfiles.tag  = 'srcfiles';
srcfiles.filter = '.*';
srcfiles.num  = [0 Inf];
srcfiles.help = {'Select files to move.'};

targetdir.type = 'files';
targetdir.name = 'Target directory';
targetdir.tag  = 'targetdir';
targetdir.filter = 'dir';
targetdir.num  = 1;
targetdir.help = {'Select target directory.'};

opts.type = 'branch';
opts.name = 'Move Files';
opts.tag  = 'movefile';
opts.val  = {srcfiles,targetdir};
opts.prog = @my_movefile;
opts.vfiles = @vfiles_movefile;
opts.help = {[...
    'This facilty allows to move files in a batch. Note that moving files ' ...
    'will not make them disappear from file selection lists. Therefore one' ...
    'has to be careful not to select the original files after they have ' ...
    'been programmed to be moved.'],'',...
    ['If image files (.*img or .*nii) are selected, corresponding hdr or mat files will ' ...
     'be moved as well, if they exist.']};

%------------------------------------------------------------------------

%------------------------------------------------------------------------
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
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function try_movefile(src,dest)
% silently try to move files
try
    movefile(src,dest);
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles_movefile(varargin)
job = varargin{1};
vf={};
for k = 1:numel(job.srcfiles)
    [p n e v] = spm_fileparts(job.srcfiles{k});
    if strncmp(e,'.img',4)||strncmp(e,'.nii',4)
        vi = strfind(e,',');
        if ~isempty(vi)
            e=e(1:vi-1);
        end;
        vf{end+1} = fullfile(job.targetdir{1},[n e v]);
        if exist(fullfile(p,[n '.mat']),'file')
            vf{end+1} = fullfile(job.targetdir{1},[n '.mat']);
        end;
        if exist(fullfile(p,[n '.hdr']))
            vf{end+1}= fullfile(job.targetdir{1},[n '.hdr']);
        end;
    else
        vf{end+1} = fullfile(job.targetdir{1},[n e v]);
    end;
end;
return;

