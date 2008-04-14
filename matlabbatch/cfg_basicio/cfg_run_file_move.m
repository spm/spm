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
% $Id: cfg_run_file_move.m 1393 2008-04-14 18:53:47Z volkmar $

rev = '$Rev: 1393 $';

if isempty(job.target{1})
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
    out.files{1} = '';
else
    out.files = {};
    for k = 1:numel(job.files)
        [p n e v] = fileparts(job.files{k});
        movefile(job.files{k}, job.target{1});
        if strcmp(e,'.img') || strcmp(e,'.nii')
            try
                movefile(fullfile(p,[n '.hdr']));
                out.files{end+1} = fullfile(job.target{1}, [n '.hdr']);
                movefile(fullfile(p,[n '.mat']));
                out.files{end+1} = fullfile(job.target{1}, [n '.mat']);
            end;
        end
        out.files{end+1} = fullfile(job.target{1}, [n e v]);
    end;
end;
