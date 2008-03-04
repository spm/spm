function out = cfg_run_file_move(job)

% Move files to another directory or delete them, if no directory is
% specified.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_run_file_move.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if isempty(job.target{1})
    for k = 1:numel(job.files)
        delete(job.files{k});
    end;
    out.files{1} = '';
else
    for k = 1:numel(job.files)
        [p n e v] = fileparts(job.files{k});
        movefile(job.files{k}, job.target{1});
        out.files{k} = fullfile(job.target{1}, [n e v]);
    end;
end;
