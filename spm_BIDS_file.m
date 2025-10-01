function [BIDS_file] = spm_BIDS_file(S)
% Function to create BIDS style filenames for saving and loading data.
% FORMAT [BIDS_file] = spm_BIDS_file(S)
%
% Input Parameters:
%     S:         Struct containing input parameters
%         category = 1xn char, describing the BIDS or derivative category.
%                   e.g. 'meg' (optional);
%         description = 1xn char, describes the filename ending. e.g.
%                   'channels' (optional)
%         type = 1xn char, describes file type, or extension. e.g. '.tsv'
%                   (default: '.mat').
%         derivative = boolean, whether or not the output is a derivative
%                   (default: true).
%         detailed = boolean, whether to include ses/task/run info in output
%                   (default: true).
%         prefix = 1xn char, describes the prefix to the filename. Useful
%                   for when spm functions have added a prefix.
%         BIDS: Struct containing BIDS information
%             directory = 1xn char providing BIDS directory.
%             sub = cell array or string, containing subject names.
%             ses = cell array or string, containing session names.
%             task = cell array or string, containing task names.
%             run = cell array or string, containing run names.
%
% Output:
%     outputs: Array of structures containing directory information
%              with fields:
%              - file: Full path and filename
%              - folder: Directory path
%              - name: Full filename
%              - ext: File extension
%              - exists: Boolean, does the file exist already?
%              - + bids fields, sub, ses, task, run
%_________________________________________________________________________
%
% Further help:
%
% spm_BIDS_file is a function that takes an input of BIDS parameters, along
% with some additional specification of the way to handle those parameters,
% and provides a single, or array, of file directories. This method may be
% used within a loop that, for example, takes a list of subject names and
% updates S.BIDS.sub, or by providing the function with the list of subject
% names directly. If insufficient information is provided the method
% assumes this is intentional and provides a limited output. This may be
% useful for accessing files which are independent of the task/run, such as
% an anatomical image in an MEG study. Or, for when only a single session
% is used, and is excluded from the file specification.
% 
% Note, when S.derivative = false, this function will not create new
% folders, to avoid breaking the organisation of existing data. Otherwise,
% it will create a new folder. The function never checks if the file
% exists.
%
% Note, this function does not enforce BIDS standards 
% (found here: https://bids-specification.readthedocs.io), see spm_BIDS for
% associated methods and checks. For example, the parameter S.prefix is 
% there to make it easier to work with files produced by SPM, rather than
% maintaining the BIDS standard. 
% 
%_________________________________________________________________________

% Nicholas Alexander
% Copyright (C) 2025 Department of Imaging Neuroscience, UCL

%==========================================================================
% - S P M   B I D S   F I L E
%==========================================================================
%-CHECK inputs
%=========================================================================
%-USERS can provide varying levels of detail, which determines the output.
%--------------------------------------------------------------------------
% Default outputs
BIDS_file = struct('file', {}, 'name', {}, 'folder', {}, 'ext', {}, ...
            'exists', {}, 'sub', {}, 'ses', {}, 'task', {}, 'run', {});

folder_only = false;
if (~isfield(S,'category') || isempty(S.category))
    disp('S.category not provided. Only returning folder output');
    folder_only = true;
end
if folder_only
    S.description = 'folder_only';
elseif (~isfield(S,'description') || isempty(S.description))
    disp('S.description is empty. Only returning folder output.')
    S.description = 'folder_only';
    folder_only = true;
end
if folder_only
    S.type = 'folder_only';
elseif (~isfield(S,'type') || isempty(S.type))
    disp('Defaulting S.type to .mat')
    S.type = '.mat';
end
if (~isfield(S,'derivative') || isempty(S.derivative))
    disp('Defaulting S.derivative to true');
    S.derivative = true;
end
if (~isfield(S,'detailed') || isempty(S.detailed))
    if ~folder_only
		S.detailed = true;
    else
		disp('Defaulting S.detailed to true');
		S.detailed = true;
    end
end
if ~isfield(S,'prefix')
    S.prefix = '';
end

% Default BIDS parameters
if ~isstruct(S.BIDS)
    error('S.BIDS must be a struct.');
end
if (~isfield(S.BIDS,'directory') || isempty(S.BIDS.directory))
    warning('No S.BIDS.directory provided. Defaulting to current directory.');
    S.BIDS.directory = cd;
end

% Convert strings to cell arrays if necessary
fields_to_convert = {'sub', 'ses', 'task', 'run'};
for f = fields_to_convert
    if isfield(S.BIDS, f{1}) && ischar(S.BIDS.(f{1}))
        S.BIDS.(f{1}) = {S.BIDS.(f{1})};
    end
end

%-ITERATE over combinations of inputs
%=========================================================================
% Prepare for iteration
sub_list = S.BIDS.sub;
if isfield(S.BIDS, 'ses') && ~isempty(S.BIDS.ses)
    ses_list = S.BIDS.ses;
else
    ses_list = {''};
end
task_list = S.BIDS.task;
run_list = S.BIDS.run;
out_idx = 1;

for sub = sub_list
    for ses = ses_list
        for task = task_list
            for run = run_list
                % Build directory
                directory = S.BIDS.directory;
                if S.derivative
                    directory = fullfile(directory, 'derivatives', S.category);
                end
                if ~isempty(sub{1})
                    directory = fullfile(directory, ['sub-', sub{1}]);
                end
                if ~isempty(ses{1})
                    directory = fullfile(directory, ['ses-', ses{1}]);
                end

                if ~folder_only
                    directory = fullfile(directory, S.category);
                end

                % Ensure directory exists.
                if ~isfolder(directory)
                    if S.derivative
                        disp(['Creating folder: ', directory]);
                        mkdir(directory);
                    else
                        warning(['Folder does not exist: ', directory]);
                    end
                end

                % Build filename
                if S.detailed
                    if isempty(ses{1})
                        filename = sprintf('%ssub-%s_task-%s_run-%s_%s%s', ...
                            S.prefix, sub{1}, task{1}, run{1}, S.description, S.type);
                    else
                        filename = sprintf('%ssub-%s_ses-%s_task-%s_run-%s_%s%s', ...
                            S.prefix, sub{1}, ses{1}, task{1}, run{1}, S.description, S.type);
                    end
                else
                    if ~isempty(sub{1})
                        filename = sprintf('%ssub-%s_%s%s', S.prefix, sub{1}, S.description, S.type);
                    else
                        filename = sprintf('%s%s%s', S.prefix, S.description, S.type);
                    end
                end

                % Add to output struct
                BIDS_file(out_idx).file = fullfile(directory,filename);
                [pathstr, name, ext] = fileparts(BIDS_file(out_idx).file);
                BIDS_file(out_idx).folder = pathstr;
                BIDS_file(out_idx).name = name;
                BIDS_file(out_idx).ext = ext;
                BIDS_file(out_idx).sub = sub;
                BIDS_file(out_idx).ses = ses;
                BIDS_file(out_idx).task = task;
                BIDS_file(out_idx).run = run;

                % Does the file exist?
                if exist(BIDS_file(out_idx).file, 'file')
                    BIDS_file(out_idx).exists = true;
                else
                    BIDS_file(out_idx).exists = false;
                end

                out_idx = out_idx + 1;
            end
        end
    end
end
end