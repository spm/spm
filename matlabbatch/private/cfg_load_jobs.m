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
% $Id: cfg_load_jobs.m 1862 2008-06-30 14:12:49Z volkmar $

rev = '$Rev: 1862 $'; %#ok

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
                cfg_message('matlabbatch:initialise:xml','LoadXML failed: ''%s''',filenames{cf});
            end;
        case '.mat'
            try
                S=load(filenames{cf});
                matlabbatch = S.matlabbatch;
            catch
                cfg_message('matlabbatch:initialise:mat','Load failed: ''%s''',filenames{cf});
            end;
        case '.m'
            opwd = pwd;
            try
                if ~isempty(p)
                    cd(p);
                end;
                clear(nam);
                eval(nam);
            catch
                cfg_message('matlabbatch:initialise:m','Eval failed: ''%s''',filenames{cf});
            end;
            cd(opwd);
            if ~exist('matlabbatch','var')
                cfg_message('matlabbatch:initialise:m','No matlabbatch job found in ''%s''', filenames{cf});
            end;
        otherwise
            cfg_message('matlabbatch:initialise:unknown','Unknown extension: ''%s''', filenames{cf});
    end;
    if exist('matlabbatch','var')
        newjobs = {newjobs{:} matlabbatch};
        clear matlabbatch;
    end;
end;
