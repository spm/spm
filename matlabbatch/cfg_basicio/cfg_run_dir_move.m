function out = cfg_run_dir_move(job)

% Move, copy or delete directory
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_dir_move.m 6920 2016-11-02 14:48:01Z guillaume $

rev = '$Rev: 6920 $'; %#ok

action = fieldnames(job.action);
action = action{1};
% delete
if strcmp(action, 'delete')
    rmdir(char(job.dir),'s');
    out = [];
else
    [dummy, dirname] = fileparts(char(job.dir));
    if strncmp(action, 'copy', 4)
        % copy
        cmd = @copyfile;
        if strcmp(action,'copyto')
            tgt = job.action.copyto{1};
        else
            tgt = job.action.(action).copyto{1};
        end
        tgt = fullfile(tgt, dirname);
        out.dir = {tgt};
    else
        % move
        cmd = @movefile;
        if strcmp(action,'moveto')
            tgt = job.action.moveto{1};
        else
            tgt = job.action.(action).moveto{1};
        end
        if exist(tgt, 'dir')
            out.dir = {fullfile(tgt, dirname)};
        else
            out.dir = tgt;
        end
    end
    
    feval(cmd, char(job.dir), tgt);
end
