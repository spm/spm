function varargout = spm_jobman(varargin)
% Main interface for SPM Batch System
% This function provides a compatibility layer between SPM and matlabbatch.
%
% FORMAT spm_jobman('initcfg')
% Initialise jobs configuration and set MATLAB path accordingly.
%
% FORMAT spm_jobman('run',job[,input1,...inputN])
% FORMAT output_list = spm_jobman('run',job[,input1,...inputN])
% FORMAT [output_list, hjob] = spm_jobman('run',job[,input1,...inputN])
% Run specified job.
% job         - filename of a job (.m or .mat), or
%               cell array of filenames, or
%               'jobs'/'matlabbatch' variable, or
%               cell array of 'jobs'/'matlabbatch' variables.
% input1,...  - optional list of input arguments. These are filled into
%               open inputs ('X->' marked in the GUI) before a job is
%               run. When using an "{:}" subscript on a cell array,
%               MATLAB expands this cell array into a comma separated 
%               list of arguments. Therefore, one can collect input
%               arguments in the right order into a cell array named e.g.
%               input_array and call spm_jobman('run',job,input_array{:})
%               to run a job using the collected inputs. For files or text
%               entry items, these inputs are the values one would specify
%               in the GUI. For menus, the item number has to be entered
%               (neither the GUI label nor the associated value that is
%               saved in a batch).
% output_list - cell array containing the output arguments from each
%               module in the job. The format and contents of these
%               outputs is defined in the configuration of each module
%               (.prog and .vout callbacks).
% hjob        - harvested job after it has been filled and run. Note that
%               this job has no dependencies any more. If one loads this
%               job back to the batch system and changes some of the
%               inputs, changed outputs will not be passed on.
%
% FORMAT job_id = spm_jobman
%        job_id = spm_jobman('interactive',job[,node])
%        job_id = spm_jobman('interactive','',node)
% Run the user interface in interactive mode.
% job         - filename of a job (.m or .mat), or
%               cell array of filenames, or
%               'jobs'/'matlabbatch' variable, or
%               cell array of 'jobs'/'matlabbatch' variables.
% node        - indicate which part of the configuration is to be used.
%               For example, it could be 'spm.spatial.coreg.estimate'.
% job_id      - can be used to manipulate this job in cfg_util. Note that
%               changes to the job in cfg_util will not show up in cfg_ui
%               unless 'Update View' is called.
%
% FORMAT jobs = spm_jobman('convert',jobs)
% Convert older batch jobs to latest version
% jobs        - char or cell array of filenames, or
%               'jobs'/'matlabbbatch' variable
%__________________________________________________________________________
% Copyright (C) 2005-2016 Wellcome Trust Centre for Neuroimaging

% Volkmar Glauche
% $Id: spm_jobman.m 7744 2019-12-03 12:38:47Z guillaume $


%__________________________________________________________________________
%
% Programmers help:
%
% FORMAT [output_list hjob] = spm_jobman('serial')
%        [output_list hjob] = spm_jobman('serial',job[,'',   input1,...inputN])
%        [output_list hjob] = spm_jobman('serial',job ,node[,input1,...inputN])
%        [output_list hjob] = spm_jobman('serial',''  ,node[,input1,...inputN])
% Run the user interface in serial mode. If job is not empty, then node
% is silently ignored. Inputs can be a list of arguments. These are passed
% on to the open inputs of the specified job/node. Each input should be
% suitable to be assigned to item.val{1}. For cfg_repeat/cfg_choice items,
% input should be a cell list of indices input{1}...input{k} into
% item.value. See cfg_util('filljob',...) for details.
%
% FORMAT spm_jobman('help',node)
%        spm_jobman('help',node,width)
% Create a cell array containing help information.  This is justified
% to be 'width' characters wide. e.g.
%     h = spm_jobman('help','spm.spatial.coreg.estimate');
%     for i=1:numel(h), fprintf('%s\n',h{i}); end
%
% FORMAT [tag, job] = spm_jobman('harvest', job_id|job|cfg_item|cfg_struct)
% Take the job with id job_id in cfg_util and extract what is
% needed to save it as a batch job (for experts only). If a (partial) job
% is given instead, the output job is augmented with default settings. 
% If the argument is a cfg_item or cfg_struct tree, it will be harvested
% outside cfg_util.  
% tag - tag of the root node of the current job/cfg_item tree
% job - harvested data from the current job/cfg_item tree
%__________________________________________________________________________
% 
% This code is based on earlier versions by John Ashburner, Philippe 
% Ciuciu and Guillaume Flandin.
% It now relies on matlabbatch
%                http://sourceforge.net/projects/matlabbatch/
% Copyright (C) 2008 Freiburg Brain Imaging
%__________________________________________________________________________


%-Force jobs configuration initialisation if needed
%--------------------------------------------------------------------------
persistent isInitCfg;
if isempty(isInitCfg) &&  ~(nargin == 1 && strcmpi(varargin{1},'initcfg'))
    % Run spm_jobman('initcfg') beforehand.
    fprintf('Initialising batch system... ');
    spm_jobman('initcfg');
    fprintf('done.\n');
end
isInitCfg = true;

%-Open GUI when called without input arguments
%--------------------------------------------------------------------------
if ~nargin
    spm_jobman('interactive');
    if nargout > 0, varargout = {findobj(0,'tag','cfg_ui')}; end
    return;
end

%-Warn about deprecated syntax
%--------------------------------------------------------------------------
action = lower(varargin{1});
switch action
    case 'run_nogui'
        warning('spm:spm_jobman:deprecated', ...
            'Callback ''%s'' is deprecated. Use ''run'' instead.',action);
        action = 'run';
    case {'spm5tospm8','spm5tospm8bulk'}
        warning('spm:spm_jobman:deprecated', ...
            'Callback ''%s'' is deprecated. Use ''convert'' instead.',action);
        action = 'convert';
end

%-Load and convert batch jobs
%--------------------------------------------------------------------------
if ismember(action, {'interactive','run','serial'})
    if nargin > 1
        % sort out job/node arguments for interactive, serial, run cmds
        if nargin>=2 && ~isempty(varargin{2})
            % do not consider node if job is given
            if ischar(varargin{2}) || iscellstr(varargin{2})
                jobs = load_jobs(varargin{2});
            elseif iscell(varargin{2})
                if iscell(varargin{2}{1})
                    % assume varargin{2} is a cell of jobs
                    jobs = varargin{2};
                else
                    % assume varargin{2} is a single job
                    jobs{1} = varargin{2};
                end
            end
            mljob = canonicalise_jobs(jobs);
        elseif ismember(action, {'interactive','serial'}) && nargin>=3 && isempty(varargin{2})
            % Node spec only allowed for 'interactive', 'serial'
            mod_cfg_id = cfg_util('tag2mod_cfg_id',varargin{3});
        else
            error('spm:spm_jobman:invalidSyntax', ...
                'Don''t know how to handle this ''%s'' call.', action);
        end
    end
end

%-Perform action
%--------------------------------------------------------------------------
switch action
    case {'initcfg'}
        if ~isdeployed
            addpath(fullfile(spm('Dir'),'matlabbatch'));
            addpath(fullfile(spm('Dir'),'config'));
        end
        try
            spm_select('init');
        catch
            S = which('spm_select','-all');
            if iscell(S) && numel(S) > 1
                fprintf('spm_select appears several times in your MATLAB path:\n');
                for i=1:numel(S)
                    if i==1
                        fprintf('  %s (SHADOWING)\n',S{1});
                    else
                        fprintf('  %s\n',S{i});
                    end
                end
            end
            rethrow(lasterror);
        end
        cfg_get_defaults('cfg_util.genscript_run', @genscript_run);
        cfg_util('initcfg'); % This must be the first call to cfg_util
        %if ~spm('cmdline')
        %    f = cfg_ui('Visible','off'); % Create invisible batch ui
        %    f0 = findobj(f, 'Tag','MenuFile'); % Add entries to file menu
        %    f2 = uimenu(f0,'Label','xxx', 'Callback',@xxx, ...
        %        'HandleVisibility','off', 'tag','xxx');
        %end
        
    case {'interactive'}
        if exist('mljob', 'var')
            cjob = cfg_util('initjob', mljob);
        elseif exist('mod_cfg_id', 'var')
            if isempty(mod_cfg_id)
                warning('spm:spm_jobman:NodeNotFound', ...
                    ['Can not find executable node ''%s'' - running '...
                    'matlabbatch without default node.'], varargin{3});
                cjob = cfg_util('initjob');
            else
                cjob = cfg_util('initjob');
                mod_job_id = cfg_util('addtojob', cjob, mod_cfg_id);
                cfg_util('harvest', cjob, mod_job_id);
            end
        else
            cjob = cfg_util('initjob');
        end
        f = findobj(0,'tag','cfg_ui');
        if isempty(f), f = cfg_ui; end
        cfg_ui('local_showjob', f, cjob);
        if nargout > 0
            varargout{1} = cjob;
        end
        
    case {'serial'}
        if exist('mljob', 'var')
            cjob = cfg_util('initjob', mljob);
        else
            cjob = cfg_util('initjob');
            if nargin > 2
                [mod_cfg_id, item_mod_id] = cfg_util('tag2cfg_id', lower(varargin{3}));
                cfg_util('addtojob', cjob, mod_cfg_id);
            end
        end
        sts = fill_run_job('serial', cjob, varargin{4:end});
        if sts
            if nargout > 0
                varargout{1} = cfg_util('getalloutputs', cjob);
            end
            if nargout > 1
                [un, varargout{2}] = cfg_util('harvestrun', cjob);
            end
            cfg_util('deljob', cjob);
        end
        
    case {'run'}
        if ~exist('mljob', 'var')
            error('Not enough input arguments.');
        end
        cjob = cfg_util('initjob', mljob);
        sts = fill_run_job('run', cjob, varargin{3:end});
        if sts
            if nargout > 0
                varargout{1} = cfg_util('getalloutputs', cjob);
            end
            if nargout > 1
                [un, varargout{2}] = cfg_util('harvestrun', cjob);
            end
            cfg_util('deljob', cjob);
        end
        
    case {'convert'}
        varargout{1} = convert_jobs(varargin{2:end});
        
    case {'harvest'}
        if nargin == 1
            error('spm:spm_jobman:CantHarvest', ...
                ['Can not harvest job without job_id. Please use ' ...
                'spm_jobman(''harvest'', job_id).']);
        elseif cfg_util('isjob_id', varargin{2})
            [tag, job] = cfg_util('harvest', varargin{2});
        elseif iscell(varargin{2})
            cjob      = cfg_util('initjob', varargin{2});
            [tag, job] =  cfg_util('harvest', cjob);
            cfg_util('deljob', cjob);
        elseif isa(varargin{2}, 'cfg_item')
            [tag, job] = harvest(varargin{2}, varargin{2}, false, false);
        elseif isstruct(varargin{2})
            % try to convert into class before harvesting
            c = cfg_struct2cfg(varargin{2});
            [tag, job] = harvest(c,c,false,false);
        else
            error('spm:spm_jobman:CantHarvestThis', ...
                'Can not harvest this argument.');
        end
        varargout{1} = tag;
        varargout{2} = job;
        
    case {'help'}
        if (nargin < 2) || isempty(varargin{2})
            node = 'spm';
        else
            node = varargin{2};
        end
        if nargin < 3
            width = 60;
        else
            width = varargin{3};
        end
        varargout{1} = cfg_util('showdocwidth', width, node);
        
    otherwise
        error('spm:spm_jobman:unknownOption','Unknown option "%s".',varargin{1});
end


%==========================================================================
% function newjobs = load_jobs(job)
%==========================================================================
function newjobs = load_jobs(job)
% Load a list of possible job files, return a cell list of jobs.
% If a job file failed to load, an empty cell is returned in the list.
filenames = cellstr(job);
newjobs = {};
for i = 1:numel(filenames)
    switch spm_file(filenames{i},'ext')
        case 'mat'
            try
                S = load(filenames{i});
                if isfield(S,'matlabbatch')
                    matlabbatch = S.matlabbatch;
                elseif isfield(S,'jobs')
                    jobs = S.jobs;
                end
            catch
                warning('spm:spm_jobman:loadFailed','Load failed: ''%s''',filenames{i});
            end
        case 'm'
            try
                str = fileread(filenames{i});
                eval(str);
            catch
                warning('spm:spm_jobman:loadFailed','Load failed: ''%s''',filenames{i});
            end
        case 'json'
            try
                S = spm_jsonread(filenames{i});
                if isstruct(S)
                    for j=1:numel(S)
                        matlabbatch{j} = S(j);
                    end
                else
                    matlabbatch = S;
                end
            catch
                warning('spm:spm_jobman:loadFailed','Load failed: ''%s''',filenames{i});
            end
        otherwise
            warning('spm:spm_jobman:unknownExt','Unknown extension: ''%s''',filenames{i});
    end
    if exist('jobs','var')
        newjobs = [newjobs(:); {jobs}];
        clear jobs;
    elseif exist('matlabbatch','var')
        newjobs = [newjobs(:); {matlabbatch}];
        clear matlabbatch;
    else
        warning('spm:spm_jobman:jobNotFound','No batch job found in ''%s''',filenames{i});
        newjobs = [newjobs(:); {[]}];
    end
end


%==========================================================================
% function varargout = convert_jobs(varargin)
%==========================================================================
function varargout = convert_jobs(varargin)
% Convert a list of jobs to latest version
if nargin && iscell(varargin{1}) && ~iscellstr(varargin{1})
    varargout = canonicalise_jobs(varargin);
    return;
elseif ~nargin || isempty(varargin{1})
    [fname, sts] = spm_select([1 Inf], 'batch', 'Select job file(s)');
    if ~sts, return; end
else
    fname = varargin{1};
end
fname     = cellstr(fname);
joblist   = load_jobs(fname);
SPMver    = spm('Ver');
outnames  = cell(numel(fname),1);
for i=1:numel(fname)
    if ~isempty(joblist{i})
        outnames{i} = spm_file(fname{i},'prefix',[lower(SPMver) '_']);
        fprintf('Initial job: %s\n', fname{i});
        fprintf('%s job: %s\n',SPMver, outnames{i});
        cjob = cfg_util('initjob', canonicalise_jobs(joblist(i)));
        cfg_util('savejob', cjob, outnames{i});
        cfg_util('deljob', cjob);
    end
end
varargout = {outnames};


%==========================================================================
% function [mljob, comp] = canonicalise_jobs(job)
%==========================================================================
function [mljob, comp] = canonicalise_jobs(job)
% job: a cell list of job data structures.
% Check whether job is a SPM5 or matlabbatch job. In the first case, all
% items in job{:} should have a fieldname of either 'temporal', 'spatial',
% 'stats', 'tools' or 'util'. If this is the case, then job will be
% assigned to mljob{1}.spm, which is the tag of the SPM root configuration
% item.
comp  = true(size(job));
mljob = cell(size(job));
for i = 1:numel(job)
    for j = 1:numel(job{i})
        comp(i) = comp(i) && any(strcmp(fieldnames(job{i}{j}), ...
            {'temporal', 'spatial', 'stats', 'tools', 'util'}));
        if ~comp(i)
            break;
        end
    end
    if comp(i)
        tmp = sub_canonicalise_job(job{i});
        for j=1:numel(tmp)
            mljob{i}{j}.spm = tmp{j};
        end
    else
        mljob{i} = job{i};
    end
end


%==========================================================================
% function njobs = sub_canonicalise_job(jobs)
%==========================================================================
function njobs = sub_canonicalise_job(jobs)
decel = struct('spatial',struct('realign',[],'coreg',[],'normalise',[]),...
               'temporal',[],...
               'stats',[],...
               'meeg',[],...
               'util',[],...
               'tools',struct('dartel',[]));
njobs  = {};
for i0 = 1:numel(jobs)
    tmp0  = fieldnames(jobs{i0});
    tmp0  = tmp0{1};
    if any(strcmp(tmp0,fieldnames(decel)))
        for i1=1:numel(jobs{i0}.(tmp0))
            tmp1  = fieldnames(jobs{i0}.(tmp0){i1});
            tmp1  = tmp1{1};
            if ~isempty(decel.(tmp0))
                if any(strcmp(tmp1,fieldnames(decel.(tmp0))))
                    for i2=1:numel(jobs{i0}.(tmp0){i1}.(tmp1))
                        njobs{end+1} = struct(tmp0,struct(tmp1,jobs{i0}.(tmp0){i1}.(tmp1){i2}));
                    end
                else
                    njobs{end+1} = struct(tmp0,jobs{i0}.(tmp0){i1});
                end
            else
                njobs{end+1} = struct(tmp0,jobs{i0}.(tmp0){i1});
            end
        end
    else
        njobs{end+1} = jobs{i0};
    end
end


%==========================================================================
% function sts = fill_run_job(action, cjob, varargin)
%==========================================================================
function sts = fill_run_job(action, cjob, varargin)
switch lower(action)
    case 'serial'
        sts = cfg_util('filljobui', cjob, @serial_ui, varargin{:});
    case 'run'
        sts = cfg_util('filljob', cjob, varargin{:});
end
if sts
    cfg_util('run', cjob);
else
    cfg_util('deljob', cjob);
    error('spm:spm_jobman:jobNotFilled', 'No executable modules, but still unresolved dependencies or incomplete module inputs.');
end


%==========================================================================
% function [val sts] = serial_ui(item)
%==========================================================================
function [val, sts] = serial_ui(item)
% wrapper function to translate cfg_util('filljobui'... input requests into
% spm_input/cfg_select calls.
sts = true;
switch class(item)
    case 'cfg_choice'
        labels = cell(size(item.values));
        values = cell(size(item.values));
        for k = 1:numel(item.values)
            labels{k} = item.values{k}.name;
            values{k} = k;
        end
        val = spm_input(item.name, 1, 'm', labels, values);
    case 'cfg_menu'
        val = spm_input(item.name, 1, 'm', item.labels, item.values);
        val = val{1};
    case 'cfg_repeat'
        labels = cell(size(item.values));
        values = cell(size(item.values));
        for k = 1:numel(item.values)
            labels{k} = item.values{k}.name;
            values{k} = k;
        end
        % enter at least item.num(1) values
        for k = 1:item.num(1)
            val(k) = spm_input(sprintf('%s(%d)', item.name, k), 1, 'm', ...
                               labels, values);
        end
        % enter more (up to varargin{3}(2) values
        labels = {labels{:} 'Done'};
        % values is a cell list of natural numbers, use -1 for Done
        values = {values{:} -1}; 
        while numel(val) < item.num(2)
            val1 = spm_input(sprintf('%s(%d)', item.name, numel(val)+1), 1, ...
                             'm', labels, values);
            if val1{1} == -1
                break;
            else
                val(end+1) = val1;
            end
        end
    case 'cfg_entry'
        val = spm_input(item.name, 1, item.strtype, '', item.num, ...
                        item.extras);
    case 'cfg_files'
        [t,sts] = cfg_getfile(item.num, item.filter, item.name, '', ...
                              item.dir, item.ufilter);
        if sts
            val = cellstr(t);
        else
            val = {};
            error('File selector was closed.');
        end
end


%==========================================================================
% function [code cont] = genscript_run
%==========================================================================
function [code, cont] = genscript_run
% Return code snippet to initialise SPM defaults and run a job generated by
% cfg_util('genscript',...) through spm_jobman.
modality = spm('CheckModality');
code{1}  = sprintf('spm(''defaults'', ''%s'');', modality);
code{2}  = 'spm_jobman(''run'', jobs, inputs{:});';
cont     = false;


%-Compatibility layer for SPM5
function varargout = interactive(varargin)
function varargout = defaults_edit(varargin)
function varargout = run_serial(varargin)
