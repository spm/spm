function newjobs = cfg_load_jobs(job)

% function newjobs = cfg_load_jobs(job)
%
% Load a list of possible job files, return a cell list of jobs. If a job
% file failed to load, an empty cell is returned in the list.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: copyright_cfg.m 218 2008-04-17 10:34:11Z glauche $

rev = '$Rev: 218 $';

if ischar(job)
    filenames = cellstr(job);
else
    filenames = job;
end;
newjobs = {};
for cf = 1:numel(filenames)
    [p,nam,ext] = fileparts(filenames{cf});
    switch ext
        case '.xml',
            try
                loadxml(filenames{cf},'matlabbatch');
            catch
                warning('cfg_util:local_load_jobs','LoadXML failed: ''%s''',filenames{cf});
            end;
        case '.mat'
            try
                S=load(filenames{cf});
                matlabbatch = S.matlabbatch;
            catch
                warning('cfg_util:local_load_jobs','Load failed: ''%s''',filenames{cf});
            end;
        case '.m'
            opwd = pwd;
            try
                if ~isempty(p)
                    cd(p);
                end;
                eval(nam);
            catch
                warning('cfg_util:local_load_jobs','Eval failed: ''%s''',filenames{cf});
            end;
            cd(opwd);
            if ~exist('matlabbatch','var')
                warning('cfg_util:local_load_jobs','No matlabbatch job found in ''%s''', filenames{cf});
            end;
        otherwise
            warning('cfg_util:local_load_jobs','Unknown extension: ''%s''', filenames{cf});
    end;
    if exist('matlabbatch','var')
        newjobs = {newjobs{:} matlabbatch};
        clear matlabbatch;
    end;
end;
