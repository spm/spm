function out = cfg_run_file_move(job)

% Move files to another directory or delete them, if no directory is
% specified. Special treatment to move .img/.hdr/.mat pairs of files
% together.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_file_move.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if isfield(job.action, 'delete')
    for k = 1:numel(job.files)
        [p n e v] = fileparts(job.files{k});
        delete(job.files{k});
        if strcmp(e,'.img') || strcmp(e,'.nii')
            try
                delete(fullfile(p,[n '.hdr']));
                delete(fullfile(p,[n '.mat']));
            end;
        end
    end;
    out = [];
else
    out.files = {};
    for k = 1:numel(job.files)
        [p n e v] = fileparts(job.files{k});
        movefile(job.files{k}, job.action.moveto{1});
        out.files{end+1} = fullfile(job.action.moveto{1}, [n e v]);
        if strcmp(e,'.img') || strcmp(e,'.nii')
            try
                movefile(fullfile(p,[n '.hdr']), job.action.moveto{1});
                out.files{end+1} = fullfile(job.action.moveto{1}, [n '.hdr']);
                movefile(fullfile(p,[n '.mat']), job.action.moveto{1});
                out.files{end+1} = fullfile(job.action.moveto{1}, [n '.mat']);
            end;
        end
    end;
end;
