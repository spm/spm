function varargout = spm_jobman(varargin)
% UI/Batching stuff
% This code is based on an earlier version by Philippe Ciuciu and
% Guillaume Flandin of Orsay, France. It provides a compatibility layer
% between SPM and matlabbatch. It translates spm_jobman callbacks into
% matlabbatch callbacks and allows to edit and run SPM5 style batch jobs.
%
% FORMAT spm_jobman
%        spm_jobman('interactive')
%        spm_jobman('interactive',job)
%        spm_jobman('interactive',job,node)
%        spm_jobman('interactive','',node)
% Runs the user interface in interactive mode.
%
% FORMAT spm_jobman('serial')
%        spm_jobman('serial',job)
%        spm_jobman('serial',job,node)
%        spm_jobman('serial','',node)
% Runs the user interface in serial mode.
%
%     node - indicates which part of the configuration is to be used.
%            For example, it could be 'jobs.spatial.coreg'.
%
%     job  - can be the name of a jobfile (as a .mat or a .xml), or a
%            'jobs' variable loaded from a jobfile.
%
% FORMAT spm_jobman('run',job)
%        spm_jobman('run_nogui',job)
% Runs a job without X11 (as long as there is no graphics output from the
% job itself). The matlabbatch system does not need graphics output to run
% a job.
%
% FORMAT spm_jobman('help',node)
%        spm_jobman('help',node,width)
% Creates a cell array containing help information.  This is justified
% to be 'width' characters wide. e.g.
%     h = spm_jobman('help','jobs.spatial.coreg.estimate');
%     for i=1:numel(h),fprintf('%s\n',h{i}); end;
%
% FORMAT spm_jobman('defaults')
% Runs the interactive defaults editor.
%
% FORMAT [tag, jobs, typ] = spm_jobman('harvest')
% Take the currently configured job in cfg_util and extract what is
% needed to save it as a batch job (for experts only).
%
% FORMAT spm_jobman('pulldown')
% Creates a pulldown 'TASKS' menu in the Graphics window.
%
% FORMAT job = spm_jobman('spm5tospm8',job)
% Takes an SPM5 job structure and returns a SPM8 compatible version.
%
% FORMAT job = spm_jobman('spm5tospm8bulk',jobfiles)
% Takes a cell string with SPM5 job filenames and saves them in SPM8
% compatible format. The new job files will be MATLAB .m files and have a
% _spm8 appended to their filename.
%
% not implemented: FORMAT spm_jobman('jobhelp')
% Creates a cell array containing help information specific for a certain
% job. Help is only printed for items where job specific help is
% present. This can be used together with spm_jobman('help') to create a
% job specific manual. This feature is available only on MATLAB R14SP2
% and higher.
%
% not implemented: FORMAT spm_jobman('chmod')
% Changes the modality for the TASKS pulldown.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Copyright (C) 2008 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: spm_jobman.m 1086 2008-01-11 15:11:20Z john $


if nargin==0
    cfg_ui;
else
    switch lower(varargin{1})
        case {'interactive','serial'}
            if nargin==1
                % uncomment to clear current job
                % cfg_util('initjob');
            elseif nargin==2 || (nargin==3 && ~isempty(varargin{2}))
                % do not consider node if job is given
                mljob = canonicalise_job(varargin{2});
            elseif nargin==3 && isempty(varargin{2})
                mod_cfg_id = cfg_util('tag2mod_cfg_id',varargin{3});
            else
                warning('spm:spm_jobman:WrongUI', ...
                        'Don''t know how to handle this ''%s'' call.', lower(varargin{1}));
            end;
            if strcmpi(varargin{1}, 'interactive')
                if exist('mljob', 'var')
                    cfg_util('initjob', mljob);
                elseif exist('mod_cfg_id', 'var')
                    if isempty(mod_cfg_id)
                        warning('spm:spm_jobman:NodeNotFound', ...
                                ['Can not find executable node ''%s'' - running '...
                                 'matlabbatch without default node.'], varargin{3});
                        cfg_util('initjob');
                    else
                        cfg_util('initjob');
                        cfg_util('addtojob', mod_cfg_id);
                    end;
                end;
                cfg_ui;
            else
                if exist('mljob', 'var')
                    cfg_serial(@serial_ui, mljob)
                else
                    if nargin > 2
                        cfg_serial(@serial_ui, lower(varargin{3}));
                    else
                        cfg_serial(@serial_ui, 'jobs');
                    end;
                end;
            end;

        case {'run','run_nogui'}
            if nargin<2
                error('Nothing to run');
            end;
            [mljob comp] = canonicalise_job(varargin{2});
            if comp
                % Run a SPM5 job - force serialisation
                cfg_util('runserial',mljob);
            else
                warning('spm:spm_jobman:runparallel', ...
                        ['%s: Trying to run job in parallel. If this does ' ...
                         'not work, run this job by calling ' ...
                         'cfg_util(''runserial'', job).'], mfilename);
                cfg_util('run',mljob);
            end;
            
        case {'spm5tospm8'}
            varargout{1} = canonicalise_job(varargin{2});

        case {'spm5tospm8bulk'}
            conv_jobs(varargin{2});
            
        case {'defaults'},
            warning('spm:spm_jobman:NotImplemented', 'Not yet implemented.');
            
        case {'pulldown'}
            pulldown;

        case {'chmod'}
            warning('spm:spm_jobman:NotImplemented', 'Callback ''%s'' not implemented.', varargin{1});

        case {'help'}
            warning('spm:spm_jobman:NotImplemented', 'Not yet implemented.');

        case {'jobhelp'}
            warning('spm:spm_jobman:NotImplemented', 'Callback ''%s'' not implemented.', varargin{1});

        case {'harvest'}
            if nargin == 1
                [varargout{1:nargout}] = cfg_util('harvest');
            else
                warning('spm:spm_jobman:CantHarvestCfg', ['Can not harvest ' ...
                                    'configuration struct argument.']);
            end;

        otherwise
            error(['"' varargin{1} '" - unknown option']);
    end;
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [mljob comp] = canonicalise_job(job)
% Check whether job is a SPM5 or matlabbatch job. In the first case, all
% items in job{:} should have a fieldname of either 'temporal', 'spatial',
% 'stats', 'tools' or 'util'. If this is the case, then job will be
% assigned to mljob{1}.spm, which is the tag of the SPM root configuration
% item.

comp = true;
for k = 1:numel(job)
    comp = comp && any(strcmp(fieldnames(job{k}),{'temporal', 'spatial', ...
                        'stats', 'tools', 'util'}));
    if ~comp
        break;
    end;
end;

if comp
    mljob{1}.spm = job;
else
    mljob = job;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function conv_jobs(varargin)
% Select a list of jobs, canonicalise each of it and save as a .m file
% using gencode.
if nargin == 0 || ~iscellstr(varargin{1})
    [fname sts] = spm_select([1 Inf], 'batch', 'Select job file(s)');
    fname = cellstr(fname);
    if ~sts return; end;
else
    fname = varargin{1};
end;

joblist = load_jobs(fname);
for k = 1:numel(fname)
    if ~isempty(joblist{k})
        [p n e v] = spm_fileparts(fname{k});
        % Save new job as *_spm8.m
        newfname = fullfile(p, sprintf('%s_spm8.m', n));
        fprintf('SPM5 job: %s\nSPM8 job: %s\n', fname{k}, newfname);
        cfg_util('initjob', canonicalise_job(joblist{k}));
        cfg_util('savejob', newfname);
    end;
end;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function load_job(varargin)
% Select a single job file, canonicalise it and display it in GUI
[fname sts] = spm_select([1 1], 'batch', 'Select job file');
if ~sts return; end;

joblist = load_jobs(fname);
if ~isempty(joblist{1})
    spm_jobman('interactive',joblist{1});
end;
return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function newjobs = load_jobs(job)
% Load a list of possible job files, return a cell list of jobs. If a job
% file failed to load, an empty cell is returned in the list.
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
            spm('Pointer','Watch');
            try
                loadxml(filenames{cf},'jobs');
            catch
                warning('LoadXML failed: ''%s''',filenames{cf});
            end;
            spm('Pointer');
        case '.mat'
            try
                S=load(filenames{cf});
                jobs = S.jobs;
            catch
                warning('Load failed: ''%s''',filenames{cf});
            end;
        case '.m'
            opwd = pwd;
            try
                cd(p);
                eval(nam);
                jobs;
            catch
                warning('Eval failed: ''%s''',filenames{cf});
            end;
            cd(opwd);
        otherwise
            warning('Unknown extension: ''%s''', filenames{cf});
    end;
    if exist('jobs','var')
        newjobs = {newjobs{:} jobs};
        clear jobs;
    else
        newjobs = {newjobs{:} {}};
    end;
end;
return;

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function pulldown
fg = spm_figure('findwin','Graphics');
if isempty(fg), return; end;
set(0,'ShowHiddenHandles','on');
delete(findobj(fg,'tag','jobs'));
set(0,'ShowHiddenHandles','off');
f0 = uimenu(fg,'Label','TASKS', ...
            'HandleVisibility','off', 'tag','jobs');
f1 = uimenu(f0,'Label','BATCH', 'Callback',@cfg_ui, ...
            'HandleVisibility','off', 'tag','jobs');
f2 = uimenu(f0,'Label','Load SPM5 job', 'Callback',@load_job, ...
            'HandleVisibility','off', 'tag','jobs', 'Separator','on');
f3 = uimenu(f0,'Label','Bulk Convert SPM5 job(s)', 'Callback',@conv_jobs, ...
            'HandleVisibility','off', 'tag','jobs');

%------------------------------------------------------------------------

%------------------------------------------------------------------------
function val = serial_ui(cmclass, cmname, varargin)
% wrapper function to translate cfg_serial input requests into
% spm_input/cfg_select calls.
switch cmclass,
    case {'cfg_choice', 'cfg_menu'},
        val = spm_input(cmname, 1, 'm', varargin{1}, varargin{2});
        val = val{1};
    case 'cfg_repeat',
        % enter at least varargin{3}[1] values
        for k = 1:varargin{3}(1)
            val(k) = spm_input(sprintf('%s(%d)', cmname, k), 1, 'm', ...
                               varargin{1}, varargin{2});
        end;
        % enter more (up to varargin{3}(2) values
        labels = {varargin{1}{:} 'Done'};
        % values is a cell list of natural numbers, use -1 for Done
        values = {varargin{2}{:} -1}; 
        while numel(val) < varargin{3}(2)
            val1 = spm_input(sprintf('%s(%d)', cmname, numel(val)+1), 1, ...
                             'm', labels, values);
            if val1{1} == -1
                break;
            else
                val(end+1) = val1;
            end;
        end;
    case 'cfg_entry',
        val = spm_input(cmname, 1, varargin{1}, '', varargin{2}, ...
                        varargin{3});
    case 'cfg_files',
        [t,sts] = cfg_getfile(varargin{1}, varargin{2}, cmname, '', ...
                              varargin{3}, varargin{4});
        if sts
            val = cellstr(t);
        else
            val = {};
        end;
end;
