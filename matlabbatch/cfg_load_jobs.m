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
% $Id: cfg_load_jobs.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

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
            try
                fid = fopen(filenames{cf},'rt');
                str = fread(fid,'*char');
                fclose(fid);
                eval(str);
            catch
                cfg_message('matlabbatch:initialise:m','Eval failed: ''%s''',filenames{cf});
            end;
            if ~exist('matlabbatch','var')
                cfg_message('matlabbatch:initialise:m','No matlabbatch job found in ''%s''', filenames{cf});
            end;
        otherwise
            cfg_message('matlabbatch:initialise:unknown','Unknown extension: ''%s''', filenames{cf});
    end;
    if exist('matlabbatch','var')
        newjobs = [newjobs(:); {matlabbatch}];
        clear matlabbatch;
    end;
end;
