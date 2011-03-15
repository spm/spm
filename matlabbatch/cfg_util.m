function varargout = cfg_util(cmd, varargin)

% This is the command line interface to the batch system. It manages the
% following structures:
% * Generic configuration structure c0. This structure will be initialised
%   to an cfg_repeat with empty .values list. Each application should
%   provide an application-specific master configuration file, which
%   describes the executable module(s) of an application and their inputs.
%   This configuration will be rooted directly under the master
%   configuration node. In this way, modules of different applications can
%   be combined with each other.
%   CAVE: the root nodes of each application must have an unique tag -
%   cfg_util will refuse to add an application which has a root tag that is
%   already used by another application.
% * Job specific configuration structure cj. This structure contains the
%   modules to be executed in a job, their input arguments and
%   dependencies between them. The layout of cj is not visible to the user.
% To address executable modules and their input values, cfg_util will
% return id(s) of unspecified type. If necessary, these id(s) should be
% stored in cell arrays in a calling application, since their internal
% format may change.
%
% The commands to manipulate these structures are described below in
% alphabetical order.
%
%  cfg_util('addapp', cfg[, def])
%
% Add an application to cfg_util. If cfg is a cfg_item, then it is used
% as initial configuration. Alternatively, if cfg is a MATLAB function,
% this function is evaluated. The return argument of this function must be
% a single variable containing the full configuration tree of the
% application to be batched.
% Optionally, a defaults configuration struct or function can be supplied.
% This function must return a single variable containing a (pseudo) job
% struct/cell array which holds defaults values for configuration items.
% These defaults should be rooted at the application's root node, not at
% the overall root node. They will be inserted by calling initialise on the
% application specific part of the configuration tree.
%
%  mod_job_id = cfg_util('addtojob', job_id, mod_cfg_id)
%
% Append module with id mod_cfg_id in the cfg tree to job with id
% job_id. Returns a mod_job_id, which can be passed on to other cfg_util
% callbacks that modify the module in the job.
%
%  [mod_job_idlist, new2old_id] = cfg_util('compactjob', job_id)
%
% Modifies the internal representation of a job by removing deleted modules
% from the job configuration tree. This will invalidate all mod_job_ids and
% generate a new mod_job_idlist.
% A translation table new2old_id is provided, where 
%  mod_job_idlist = old_mod_job_idlist{new2old_id}
% translates between an old id list and the compact new id list.
%
%  cfg_util('delfromjob', job_id, mod_job_id)
%
% Delete a module from a job.
%
%  cfg_util('deljob', job_id)
%
% Delete job with job_id from the job list. 
%
%  sts = cfg_util('filljob', job_id, input1, ..., inputN)
%  sts = cfg_util('filljobui', job_id, ui_fcn, input1, ..., inputN)
%
% Fill missing inputs in a job from a list of input items. If an can not be
% filled by the specified input, this input will be discarded. If
% cfg_util('filljobui'...) is called, [val sts] = ui_fcn(item) will be run
% and should return a value which is suitable for setval(item, val,
% false). sts should be set to true if input should continue with the
% next item. This can result in an partially filled job. If ui_fcn is
% interrupted, the job will stay unfilled.
% If cfg_util('filljob'...) is called, the current job can become partially
% filled.
% Returns the all_set status of the filled job, returns always false if
% ui_fcn is interrupted.
%
%  cfg_util('gencode', fname, apptag|cfg_id[, tropts])
%
% Generate code from default configuration structure, suitable for
% recreating the tree structure. Note that function handles may not be
% saved properly. By default, the entire tree is saved into a file fname.
% If tropts is given as a traversal option specification, code generation
% will be split at the nodes matching tropts.stopspec. Each of these nodes will
% generate code in a new file with filename <fname>_<tag of node>, and the
% nodes up to tropts.stopspec will be saved into fname.
% If a file named <fname>_mlb_preamble.m exists in the folder where the
% configuration code is being written, it will be read in literally
% and its contents will be prepended to each of the created files. This
% allows to automatically include e.g. copyright or revision.
%
%  cfg_util('genscript', job_id, scriptdir, filename)
% 
% Generate a script which collects missing inputs of a batch job and runs
% the job using cfg_util('filljob', ...). The script will be written to
% file filename.m in scriptdir, the job will be saved to filename_job.m in
% the same folder. The script must be completed by adding code to collect
% the appropriate inputs for the job.
%
%  outputs = cfg_util('getAllOutputs', job_id)
% 
% outputs - cell array with module outputs. If a module has not yet been
%           run, a cfg_inv_out object is returned.
%
%  voutputs = cfg_util('getAllVOutputs', job_id[, mod_job_id])
% 
% voutputs - cell array with virtual output descriptions (cfg_dep objects).
%            These describe the structure of the job outputs. To create
%            dependencies, they can be entered into matching input objects
%            in subsequent modules of the same job.
%            If mod_job_id is supplied, only virtual output descriptions of
%            the referenced module are returned.
%
%  [tag, val] = cfg_util('harvest', job_id[, mod_job_id])
%
% Harvest is a method defined for all 'cfg_item' objects. It collects the
% entered values and dependencies of the input items in the tree and
% assembles them in a struct/cell array.
% If no mod_job_id is supplied, the internal configuration tree will be
% cleaned up before harvesting. Dependencies will not be resolved in this
% case. The internal state of cfg_util is not modified in this case. The
% structure returned in val may be saved to disk as a job and can be loaded
% back into cfg_util using the 'initjob' command.
% If a mod_job_id is supplied, only the relevant part of the configuration
% tree is harvested, dependencies are resolved and the internal state of
% cfg_util is updated. In this case, the val output is only part of a job
% description and can not be loaded back into cfg_util.
%
%  [tag, appdef] = cfg_util('harvestdef'[, apptag|cfg_id])
%
% Harvest the defaults branches of the current configuration tree. If
% apptag is supplied, only the subtree of that application whose root tag
% matches apptag/whose id matches cfg_id is harvested. In this case,
% appdef is a struct/cell array that can be supplied as a second argument
% in application initialisation by cfg_util('addapp', appcfg,
% appdef). 
% If no application is specified, defaults of all applications will be
% returned in one struct/cell array. 
% 
%  [tag, val] = cfg_util('harvestrun', job_id)
%
% Harvest data of a job that has been (maybe partially) run, resolving
% all dependencies that can be resolved. This can be used to document
% what has actually been done in a job and which inputs were passed to
% modules with dependencies.
% If the job has not been run yet, tag and val will be empty.
%
%  cfg_util('initcfg')
%
% Initialise cfg_util configuration. All currently added applications and
% jobs will be cleared.
% Initial application data will be initialised to a combination of
% cfg_mlbatch_appcfg.m files in their order found on the MATLAB path. Each
% of these config files should be a function with calling syntax
%   function [cfg, def] = cfg_mlbatch_appcfg(varargin) 
% This function should do application initialisation (e.g. add
% paths). cfg and def should be configuration and defaults data
% structures or the name of m-files on the MATLAB path containing these
% structures. If no defaults are provided, the second output argument
% should be empty.
% cfg_mlbatch_appcfg files are executed in the order they are found on
% the MATLAB path with the one first found taking precedence over
% following ones.
%
%  cfg_util('initdef', apptag|cfg_id[, defvar])
%
% Set default values for application specified by apptag or
% cfg_id. If defvar is supplied, it should be any representation of a
% defaults job as returned by cfg_util('harvestdef', apptag|cfg_id),
% i.e. a MATLAB variable, a function creating this variable...
% Defaults from defvar are overridden by defaults specified in .def
% fields.
% New defaults only apply to modules added to a job after the defaults
% have been loaded. Saved jobs and modules already present in the current
% job will not be changed.
%
%  [job_id, mod_job_idlist] = cfg_util('initjob'[, job])
%
% Initialise a new job. If no further input arguments are provided, a new
% job without modules will be created.
% If job is given as input argument, the job tree structure will be
% loaded with data from the struct/cell array job and a cell list of job
% ids will be returned.
% The new job will be appended to an internal list of jobs. It must
% always be referenced by its job_id.
%
%  sts = cfg_util('isjob_id', job_id)
%  sts = cfg_util('ismod_cfg_id', mod_cfg_id)
%  sts = cfg_util('ismod_job_id', job_id, mod_job_id)
%  sts = cfg_util('isitem_mod_id', item_mod_id)
% Test whether the supplied id seems to be of the queried type. Returns
% true if the id matches the data format of the queried id type, false
% otherwise. For item_mod_ids, no checks are performed whether the id is
% really valid (i.e. points to an item in the configuration
% structure). This can be used to decide whether 'list*' or 'tag2*'
% callbacks returned valid ids.
%
%  [mod_cfg_idlist, stop, [contents]] = cfg_util('listcfg[all]', mod_cfg_id, find_spec[, fieldnames])
%
% List modules and retrieve their contents in the cfg tree, starting at
% mod_cfg_id. If mod_cfg_id is empty, search will start at the root level
% of the tree. The returned mod_cfg_id_list is always relative to the root
% level of the tree, not to the mod_cfg_id of the start item. This search
% is designed to stop at cfg_exbranch level. Its behaviour is undefined if
% mod_cfg_id points to an item within an cfg_exbranch. See 'match' and
% 'cfg_item/find' for details how to specify find_spec. A cell list of
% matching modules is returned.
% If the 'all' version of this command is used, also matching
% non-cfg_exbranch items up to the first cfg_exbranch are returned. This
% can be used to build a menu system to manipulate configuration.
% If a cell array of fieldnames is given, contents of the specified fields
% will be returned. See 'cfg_item/list' for details. This callback is not
% very specific in its search scope. To find a cfg_item based on the
% sequence of tags of its parent items, use cfg_util('tag2mod_cfg_id',
% tagstring) instead.
%
%  [item_mod_idlist, stop, [contents]] = cfg_util('listmod', job_id, mod_job_id, item_mod_id, find_spec[, tropts][, fieldnames])
%  [item_mod_idlist, stop, [contents]] = cfg_util('listmod', mod_cfg_id, item_mod_id, find_spec[, tropts][, fieldnames])
%
% Find configuration items starting in module mod_job_id in the job
% referenced by job_id or in module mod_cfg_id in the defaults tree,
% starting at item item_mod_id. If item_mod_id is an empty array, start
% at the root of a module. By default, search scope are the filled items
% of a module. See 'match' and 'cfg_item/find' for details how to specify
% find_spec and tropts and how to search the default items instead of the
% filled ones. A cell list of matching items is returned.
% If a cell array of fieldnames is given, contents of the specified fields
% will be returned. See 'cfg_item/list' for details.
%
%  sts = cfg_util('match', job_id, mod_job_id, item_mod_id, find_spec)
%
% Returns true if the specified item matches the given find spec and false
% otherwise. An empty item_mod_id means that the module node itself should
% be matched.
%
%  new_mod_job_id = cfg_util('replicate', job_id, mod_job_id[, item_mod_id, val])
%
% If no item_mod_id is given, replicate a module by appending it to the
% end of the job with id job_id. The values of all items will be
% copied. This is in contrast to 'addtojob', where a module is added with
% default settings. Dependencies where this module is a target will be
% kept, whereas source dependencies will be dropped from the copied module.
% If item_mod_id points to a cfg_repeat object within a module, its
% setval method is called with val. To achieve replication, val(1) must
% be finite and negative, and val(2) must be the index into item.val that
% should be replicated. All values are copied to the replicated entry.
%
%  cfg_util('run'[, job|job_id])
%
% Run the currently configured job. If job is supplied as argument and is
% a harvested job, then cfg_util('initjob', job) will be called first. If
% job_id is supplied and is a valid job_id, the job with this job id will
% be run.
% The job is harvested and dependencies are resolved if possible.
% If cfg_get_defaults('cfg_util.runparallel') returns true, all
% modules without unresolved dependencies will be run in arbitrary order.
% Then the remaining modules are harvested again and run, if their
% dependencies can be resolved. This process is iterated until no modules
% are left or no more dependencies can resolved. In a future release,
% independent modules may run in parallel, if there are licenses to the
% Distributed Computing Toolbox available.
% Note that this requires dependencies between modules to be described by
% cfg_dep objects. If a module e.g. relies on file output of another module
% and this output is already specified as a filename of a non-existent
% file, then the dependent module may be run before the file is created.
% Side effects (changes in global variables, working directories) are
% currently not modeled by dependencies.
% If a module fails to execute, computation will continue on modules that
% do not depend on this module. An error message will be logged and the
% module will be reported as 'failed to run' in the MATLAB command window.
%
%  cfg_util('runserial'[, job|job_id])
%
% Like 'run', but force cfg_util to run the job as if each module was
% dependent on its predecessor. If cfg_get_defaults('cfg_util.runparallel')
% returns false, cfg_util('run',...) and cfg_util('runserial',...) are
% identical.
%
%  cfg_util('savejob', job_id, filename)
%
% The current job will be save to the .m file specified by filename. This
% .m file contains MATLAB script code to recreate the job variable. It is
% based on gencode (part of this MATLAB batch system) for all standard
% MATLAB types. For objects to be supported, they must implement their own
% gencode method.
%
%  cfg_util('savejobrun', job_id, filename)
%
% Save job after it has been run, resolving dependencies (see
% cfg_util('harvestrun',...)). If the job has not been run yet, nothing
% will be saved.
%
%  sts = cfg_util('setval', job_id, mod_job_id, item_mod_id, val)
%
% Set the value of item item_mod_id in module mod_job_id to val. If item is
% a cfg_choice, cfg_repeat or cfg_menu and val is numeric, the value will
% be set to item.values{val(1)}. If item is a cfg_repeat and val is a
% 2-vector, then the min(val(2),numel(item.val)+1)-th value will be set
% (i.e. a repeat added or replaced). If val is an empty cell, the value of
% item will be cleared.
% sts returns the status of all_set_item after the value has been
% set. This can be used to check whether the item has been successfully
% set.
% Once editing of a module has finished, the module needs to be harvested
% in order to update dependencies from and to other modules.
%
%  cfg_util('setdef', mod_cfg_id, item_mod_id, val)
% 
% Like cfg_util('setval',...) but set items in the defaults tree. This is
% only supported for cfg_leaf items, not for cfg_choice, cfg_repeat,
% cfg_branch items.
% Defaults only apply to new jobs, not to already configured ones.
%
%  doc = cfg_util('showdoc', tagstr|cfg_id|(job_id, mod_job_id[, item_mod_id]))
%
% Return help text for specified item. Item can be either a tag string or
% a cfg_id in the default configuration tree, or a combination of job_id,
% mod_job_id and item_mod_id from the current job.
% The text returned will be a cell array of strings, each string
% containing one paragraph of the help text. In addition to the help
% text, hints about valid values, defaults etc. are displayed.
%
%  doc = cfg_util('showdocwidth', handle|width, tagstr|cfg_id|(job_id, mod_job_id[, item_mod_id]))
%
% Same as cfg_util('showdoc', but use handle or width to determine the
% width of the returned strings.
%
%  [mod_job_idlist, str, sts, dep sout] = cfg_util('showjob', job_id[, mod_job_idlist])
%
% Return information about the current job (or the part referenced by the
% input cell array mod_job_idlist). Output arguments
% * mod_job_idlist - cell list of module ids (same as input, if provided)
% * str            - cell string of names of modules 
% * sts            - array of all set status of modules
% * dep            - array of dependency status of modules
% * sout           - array of output description structures 
% Each module configuration may provide a callback function 'vout' that
% returns a struct describing module output variables. See 'cfg_exbranch'
% for details about this callback, output description and output structure.
% The module needs to be harvested before to make output_struct available.
% This information can be used by the calling application to construct a
% dependency object which can be passed as input to other modules. See
% 'cfg_dep' for details about dependency objects.
%
%  [mod_cfg_id, item_mod_id] = cfg_util('tag2cfg_id', tagstr)
%
% Return a mod_cfg_id for the cfg_exbranch item that is the parent to the
% item in the configuration tree whose parents have tag names as in the
% dot-delimited tag string. item_mod_id is relative to the cfg_exbranch
% parent. If tag string matches a node above cfg_exbranch level, then
% item_mod_id will be invalid and mod_cfg_id will point to the specified
% node.
% Use cfg_util('ismod_cfg_id') and cfg_util('isitem_mod_id') to determine
% whether returned ids are valid or not.
% Tag strings should begin at the root level of an application configuration, 
% not at the matlabbatch root level.
%
%  mod_cfg_id = cfg_util('tag2mod_cfg_id', tagstr)
%
% Same as cfg_util('tag2cfg_id', tagstr), but it only returns a proper
% mod_cfg_id. If none of the tags in tagstr point to a cfg_exbranch, then
% mod_cfg_id will be invalid.
%
% The layout of the configuration tree and the types of configuration items
% have been kept compatible to a configuration system and job manager
% implementation in SPM5 (Statistical Parametric Mapping, Copyright (C)
% 2005 Wellcome Department of Imaging Neuroscience). This code has been
% completely rewritten based on an object oriented model of the
% configuration tree.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_util.m 4252 2011-03-15 17:49:51Z volkmar $

rev = '$Rev: 4252 $';

%% Initialisation of cfg variables
% load persistent configuration data, initialise if necessary

% generic configuration structure
persistent c0;
% job specific configuration structure
% This will be initialised to a struct (array) with fields cj and
% id2subs. When initialising a new job, it will be appended to this
% array. Jobs in this array may be cleared by setting cj and id2subs to
% [].
% field cj:
% configuration tree of this job.
% field cjid2subs:
% cell array that maps ids to substructs into the configuration tree -
% ids do not change for a cfg_util life time, while the actual position
% of a module in cj may change due to adding/removing modules. This would
% also allow to reorder modules in cj without changing their id.
persistent jobs;

if isempty(c0) && ~strcmp(cmd,'initcfg')
    % init, if not yet done
    cfg_util('initcfg');
end

%% Callback switches
% evaluate commands
switch lower(cmd),
    case 'addapp',
        [c0 jobs] = local_addapp(c0, jobs, varargin{:});
    case 'addtojob',
        cjob = varargin{1};
        mod_cfg_id = varargin{2};
        if cfg_util('isjob_id', cjob) && cfg_util('ismod_cfg_id', mod_cfg_id)
            [jobs(cjob), mod_job_id] = local_addtojob(jobs(cjob), mod_cfg_id);
            varargout{1} = mod_job_id;
        end
    case 'compactjob',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob)
            [jobs(cjob), n2oid] = local_compactjob(jobs(cjob));
            varargout{1} = num2cell(1:numel(jobs(cjob).cjid2subs));
            varargout{2} = n2oid;
        end
    case 'delfromjob',
        cjob = varargin{1};
        mod_job_id = varargin{2};
        if cfg_util('ismod_job_id', cjob, mod_job_id)
            jobs(cjob) = local_delfromjob(jobs(cjob), mod_job_id);
        end
    case 'deljob',
        if cfg_util('isjob_id', varargin{1})
            if varargin{1} == numel(jobs) && varargin{1} > 1
                jobs = jobs(1:end-1);
            else
                jobs(varargin{1}).cj = c0;
                jobs(varargin{1}).cjid2subs = {};
                jobs(varargin{1}).cjrun = [];
                jobs(varargin{1}).cjid2subsrun = {};
            end
        end
    case 'dumpcfg'
        %-Locate batch configs and copy them
        apps = which('cfg_mlbatch_appcfg','-all');
        appcfgs = cell(size(apps));
        p = fileparts(mfilename('fullpath'));
        for k = 1:numel(apps)
            appcfgs{k}   = fullfile(p,'private',sprintf('cfg_mlbatch_appcfg_%d.m',k));
            copyfile(apps{k}, appcfgs{k});
        end
        cmaster = fullfile(p, 'private', 'cfg_mlbatch_appcfg_master.m');
        fid     = fopen(cmaster,'w');
        fprintf(fid,'function cfg_mlbatch_appcfg_master\n');
        for k = 1:numel(apps)
            fprintf(fid,'[cfg, def] = cfg_mlbatch_appcfg_%d;\n', k);
            fprintf(fid,'cfg_util(''addapp'', cfg, def);\n');
        end
        fclose(fid);
    case 'filljob',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob)
            try
                jobs(cjob).cj = fillvals(jobs(cjob).cj, varargin(2:end), []);
                sts = all_set(jobs(cjob).cj);
            catch
                sts = false;
            end
        else
            sts = false;
        end
        varargout{1} = sts;    
    case 'filljobui',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob)
            try
                jobs(cjob).cj = fillvals(jobs(cjob).cj, varargin(3:end), varargin{2});
                sts = all_set(jobs(cjob).cj);
            catch
                sts = false;
            end
        else
            sts = false;
        end
        varargout{1} = sts;
    case 'gencode',
        fname = varargin{1};
        cm = local_getcm(c0, varargin{2});
        if nargin > 3
            tropts = varargin{3};
        else
            % default for SPM5 menu structure
            tropts = cfg_tropts(cfg_findspec, 1, 2, 0, Inf, true);
        end
        local_gencode(cm, fname, tropts);
    case 'genscript',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob)
            % Get information about job
            [mod_job_idlist, str, sts, dep sout] = cfg_util('showjob', cjob);
            opensel = find(~sts);
            fspec = cfg_findspec({{'hidden',false}});
            tropts = cfg_tropts({{'hidden','true'}},1,inf,1,inf,false);
            oclass = cell(1,numel(opensel));
            onames = cell(1,numel(opensel));
            % List modules with open inputs
            for cmind = 1:numel(opensel)
                [item_mod_idlist, stop, contents] = cfg_util('listmod', cjob, mod_job_idlist{opensel(cmind)}, [], fspec, tropts, {'class','name','all_set_item'});
                % module name is 1st in list
                % collect classes and names of open inputs
                oclass{cmind} = contents{1}(~[contents{3}{:}]);
                onames{cmind} = cellfun(@(iname)sprintf([contents{2}{1} ': %s'],iname),contents{2}(~[contents{3}{:}]),'UniformOutput',false);
            end
            % Generate filenames, save job
            scriptdir = char(varargin{2});
            [un filename]  = fileparts(varargin{3});
            scriptfile = fullfile(scriptdir, [filename '.m']);
            jobfile    = {fullfile(scriptdir, [filename '_job.m'])};
            cfg_util('savejob', cjob, char(jobfile));
            % Prepare script code
            oclass = [oclass{:}];
            onames = [onames{:}];
            % Document open inputs
            script = {'% List of open inputs'};
            for cmind = 1:numel(oclass)
                script{end+1} = sprintf('%% %s - %s', onames{cmind}, oclass{cmind});
            end
            % Create script stub code
            script{end+1} = 'nrun = X; % enter the number of runs here';
            script        = [script{:} gencode(jobfile)];
            script{end+1} = 'jobs = repmat(jobfile, 1, nrun);';
            script{end+1} = sprintf('inputs = cell(%d, nrun);', numel(oclass));
            script{end+1} = 'for crun = 1:nrun';
            for cmind = 1:numel(oclass)
                script{end+1} = sprintf('    inputs{%d, crun} = MATLAB_CODE_TO_FILL_INPUT; %% %s - %s', cmind, onames{cmind}, oclass{cmind});
            end
            script{end+1} = 'end';
            genscript_run = cfg_get_defaults('cfg_util.genscript_run');
            if ~isempty(genscript_run) && subsasgn_check_funhandle(genscript_run)
                [s1 cont] = feval(genscript_run);
                script = [script(:); s1(:)];
            else
                cont = true;
            end
            if cont
                script{end+1} = 'job_id = cfg_util(''initjob'', jobs);';
                script{end+1} = 'sts    = cfg_util(''filljob'', job_id, inputs{:});';
                script{end+1} = 'if sts';
                script{end+1} = '    cfg_util(''run'', job_id);';
                script{end+1} = 'end';
                script{end+1} = 'cfg_util(''deljob'', job_id);';
            end
            fid = fopen(scriptfile, 'wt');
            fprintf(fid, '%s\n', script{:});
            fclose(fid);
        end
    case 'getalloutputs',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob) && ~isempty(jobs(cjob).cjrun)
            varargout{1} = cellfun(@(cid)subsref(jobs(cjob).cjrun, [cid substruct('.','jout')]), jobs(cjob).cjid2subsrun, 'UniformOutput',false);
        else
            varargout{1} = {};
        end
    case 'getallvoutputs',
        if nargin == 2 && cfg_util('isjob_id', varargin{1})
            cjob = varargin{1};
            vmods = ~cellfun('isempty',jobs(cjob).cjid2subs); % Ignore deleted modules
            varargout{1} = cellfun(@(cid)subsref(jobs(cjob).cj, [cid substruct('.','sout')]), jobs(cjob).cjid2subs(vmods), 'UniformOutput',false);
        elseif nargin > 2 && cfg_util('ismod_job_id',varargin{1},varargin{2})
            cjob = varargin{1};
            cmod = varargin{2};
            varargout{1} = {subsref(jobs(cjob).cj, [jobs(cjob).cjid2subs{cmod} substruct('.','sout')])};
        else
            varargout{1} = {};
        end
    case 'harvest',
        tag = '';
        val = [];
        cjob = varargin{1};
        if nargin == 2
            if cfg_util('isjob_id', cjob)
                % harvest entire job
                % do not resolve dependencies
                cj1 = local_compactjob(jobs(cjob));
                [tag val] = harvest(cj1.cj, cj1.cj, false, false);
            end
        else
            mod_job_id = varargin{2};
            if cfg_util('ismod_job_id', cjob, mod_job_id)
                [tag val u3 u4 u5 jobs(cjob).cj] = harvest(subsref(jobs(cjob).cj, ...
                                                                  jobs(cjob).cjid2subs{mod_job_id}), ...
                                                           jobs(cjob).cj, ...
                                                           false, true);
            end
        end
        varargout{1} = tag;
        varargout{2} = val;
    case 'harvestdef',
        if nargin == 1
            % harvest all applications
            cm = c0;
        else
            cm = local_getcm(c0, varargin{1});
        end
        [tag defval] = harvest(cm, cm, true, false);
        varargout{1} = tag;
        varargout{2} = defval;
    case 'harvestrun',
        tag = '';
        val = [];
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob) && ~isempty(jobs(cjob).cjrun)
            [tag val] = harvest(jobs(cjob).cjrun, jobs(cjob).cjrun, false, ...
                                true);
        end            
        varargout{1} = tag;
        varargout{2} = val;
    case 'initcfg',
        [c0 jobs cjob] = local_initcfg;
        local_initapps;
    case 'initdef',
        [cm id] = local_getcm(c0, varargin{1});
        cm = local_initdef(cm, varargin{2});
        c0 = subsasgn(c0, id{1}, cm);
    case 'initjob'
        if isempty(jobs(end).cjid2subs)
            cjob = numel(jobs);
        else
            cjob = numel(jobs)+1;
        end
        if nargin == 1
            jobs(cjob).c0 = c0;
            jobs(cjob).cj = jobs(cjob).c0;
            jobs(cjob).cjid2subs = {};
            jobs(cjob).cjrun = [];
            jobs(cjob).cjid2subsrun = {};
            varargout{1} = cjob;
            varargout{2} = {};
            return;
        elseif ischar(varargin{1}) || iscellstr(varargin{1})
            [job jobdedup] = cfg_load_jobs(varargin{1});
        elseif iscell(varargin{1}) && iscell(varargin{1}{1})
            % try to initialise cell array of jobs
            job = varargin{1};
            jobdedup = NaN; % Postpone duplicate detection
        else
            % try to initialise single job
            job{1} = varargin{1};
            jobdedup = 1;
        end
        % job should be a cell array of job structures
        isjob = true(size(job));
        for k = 1:numel(job)
            isjob(k) = iscell(job{k}) && all(cellfun('isclass', job{k}, 'struct'));
        end
        job = job(isjob);
        if isempty(job)
            cfg_message('matlabbatch:initialise:invalid','No valid job.');
        else
            if any(isnan(jobdedup))
                % build up list of unique jobs
                jobdedup = NaN*ones(size(job));
                cu = 0;
                for k = 1:numel(job)
                    if isnan(jobdedup(k))
                        % found new candidate
                        cu = cu+1;
                        jobdedup(k) = cu;
                        % look for similar jobs under remaining candidates
                        csel = find(isnan(jobdedup));
                        eqind = cellfun(@(cjob)isequalwithequalnans(cjob,job{k}),job(csel));
                        jobdedup(csel(eqind)) = cu;
                    end
                end
            else
                jobdedup = jobdedup(isjob);
            end
            jobs(cjob).c0 = c0;
            [jobs(cjob) mod_job_idlist] = local_initjob(jobs(cjob), job, jobdedup);
            varargout{1} = cjob;
            varargout{2} = mod_job_idlist;
        end
    case 'isitem_mod_id'
        varargout{1} = isempty(varargin{1}) || ...
            (isstruct(varargin{1}) && ...
            all(isfield(varargin{1}, {'type','subs'})));
    case 'isjob_id'
        varargout{1} = ~isempty(varargin{1}) && ...
            isnumeric(varargin{1}) && ...
            varargin{1} <= numel(jobs) ...
            && (~isempty(jobs(varargin{1}).cjid2subs) ...
                || varargin{1} == numel(jobs));            
    case 'ismod_cfg_id'
        varargout{1} = isstruct(varargin{1}) && ...
            all(isfield(varargin{1}, {'type','subs'}));
    case 'ismod_job_id'
        varargout{1} = cfg_util('isjob_id', varargin{1}) && ...
            isnumeric(varargin{2}) && ...
            varargin{2} <= numel(jobs(varargin{1}).cjid2subs) ...
            && ~isempty(jobs(varargin{1}).cjid2subs{varargin{2}});
    case {'listcfg','listcfgall'}
        % could deal with hidden/modality fields here
        if strcmpi(cmd(end-2:end), 'all')
            exspec = cfg_findspec({});
        else
            exspec = cfg_findspec({{'class','cfg_exbranch'}});
        end
        % Stop traversal at hidden flag
        % If user input find_spec contains {'hidden',false}, then a hidden
        % node will not match and will not be listed. If a hidden node
        % matches, it will return with a stop-flag set.
        tropts = cfg_tropts({{'class','cfg_exbranch','hidden',true}}, 1, Inf, 0, Inf, true);
        % Find start node
        if isempty(varargin{1})
            cs = c0;
            sid = [];
        else
            cs = subsref(c0, varargin{1});
            sid = varargin{1};
        end
        if nargin < 4
            [id stop] = list(cs, [varargin{2} exspec], tropts);
            for k=1:numel(id)
                id{k} = [sid id{k}];
            end
            varargout{1} = id;
            varargout{2} = stop;
        else
            [id stop val] = list(cs, [varargin{2} exspec], tropts, varargin{3});
            for k=1:numel(id)
                id{k} = [sid id{k}];
            end
            varargout{1} = id;
            varargout{2} = stop;
            varargout{3} = val;
        end
    case 'listmod'
        if cfg_util('ismod_job_id', varargin{1}, varargin{2}) && ...
                cfg_util('isitem_mod_id', varargin{3})
            cjob        = varargin{1};
            mod_job_id  = varargin{2};
            item_mod_id = varargin{3};
            nids        = 3;
            if isempty(item_mod_id)
                cm = subsref(jobs(cjob).cj, jobs(cjob).cjid2subs{mod_job_id});
            else
                cm = subsref(jobs(cjob).cj, [jobs(cjob).cjid2subs{mod_job_id} item_mod_id]);
            end
        elseif cfg_util('ismod_cfg_id', varargin{1}) && ...
                cfg_util('isitem_mod_id', varargin{2})
            mod_cfg_id  = varargin{1};
            item_mod_id = varargin{2};
            nids        = 2;
            if isempty(varargin{2})
                cm = subsref(c0, mod_cfg_id);
            else
                cm = subsref(c0, [mod_cfg_id item_mod_id]);
            end
        end
        findspec = varargin{nids+1};
        if (nargin > nids+2 && isstruct(varargin{nids+2})) || nargin > nids+3
            tropts = varargin{nids+2};
        else
            tropts = cfg_tropts({{'hidden',true}}, 1, Inf, 0, Inf, false);
        end
        if (nargin > nids+2 && iscellstr(varargin{nids+2}))
            fn = varargin{nids+2};
        elseif nargin > nids+3
            fn = varargin{nids+3};
        else
            fn = {};
        end
        if isempty(fn)
            [id stop] = list(cm, findspec, tropts);
            varargout{1} = id;
            varargout{2} = stop;
        else
            [id stop val] = list(cm, findspec, tropts, fn);
            varargout{1} = id;
            varargout{2} = stop;
            varargout{3} = val;
        end
    case 'match'
        res = {};
        cjob        = varargin{1};
        mod_job_id  = varargin{2};
        item_mod_id = varargin{3};
        if cfg_util('ismod_job_id', cjob, mod_job_id) && ...
                cfg_util('isitem_mod_id', item_mod_id)
            if isempty(item_mod_id)
                cm = subsref(jobs(cjob).cj, jobs(cjob).cjid2subs{mod_job_id});
            else
                cm = subsref(jobs(cjob).cj, [jobs(cjob).cjid2subs{mod_job_id} item_mod_id]);
            end
            res = match(cm, varargin{4});
        end
        varargout{1} = res;
    case 'replicate'
        cjob       = varargin{1};
        mod_job_id = varargin{2};
        if cfg_util('ismod_job_id', cjob, mod_job_id)
            if nargin == 3
                % replicate module
                [jobs(cjob) id] = local_replmod(jobs(cjob), mod_job_id);
            elseif nargin == 5 && ~isempty(varargin{3}) && ...
                    cfg_util('isitem_mod_id', varargin{3})
                % replicate val entry of cfg_repeat, use setval with sanity
                % check
                cm = subsref(jobs(cjob).cj, [jobs(cjob).cjid2subs{mod_job_id}, varargin{3}]);
                if isa(cm, 'cfg_repeat')
                    cm = setval(cm, varargin{4}, false);
                    jobs(cjob).cj = subsasgn(jobs(cjob).cj, ...
                                         [jobs(cjob).cjid2subs{mod_job_id}, ...
                                        varargin{3}], cm);
                end
            end
            % clear run configuration
            jobs(cjob).cjrun = [];
            jobs(cjob).cjid2subsrun = {};
        end
    case {'run','runserial','cont','contserial'}
        if cfg_util('isjob_id',varargin{1})
            cjob = varargin{1};
            dflag = false;
        else
            cjob = cfg_util('initjob',varargin{1});
            dflag = true;
        end
        pflag = any(strcmpi(cmd, {'run','cont'})) && cfg_get_defaults([mfilename '.runparallel']);
        cflag = any(strcmpi(cmd, {'cont','contserial'}));
        [jobs(cjob) err] = local_runcj(jobs(cjob), cjob, pflag, cflag);
        if dflag
            cfg_util('deljob', cjob);
        end
        if ~isempty(err)
            cfg_message(err);
        end
    case {'savejob','savejobrun'}
        cjob = varargin{1};
        if strcmpi(cmd,'savejob')
            [tag matlabbatch] = cfg_util('harvest', cjob);
        else
            [tag matlabbatch] = cfg_util('harvestrun', cjob);
        end
        if isempty(tag)
            cfg_message('matlabbatch:cfg_util:savejob:nojob', ...
                    'Nothing to save for job #%d', cjob);
        else
            [p n e] = fileparts(varargin{2});
            switch lower(e)
                case '.mat',
                    save(varargin{2},'matlabbatch','-v6');
                case '.m'
                    jobstr = gencode(matlabbatch, tag);
                    fid = fopen(fullfile(p, [n '.m']), 'wt');
                    fprintf(fid, '%%-----------------------------------------------------------------------\n');
                    fprintf(fid, '%% Job configuration created by %s (rev %s)\n', mfilename, rev);
                    fprintf(fid, '%%-----------------------------------------------------------------------\n');
                    fprintf(fid, '%s\n', jobstr{:});
                    fclose(fid);
                otherwise
                    cfg_message('matlabbatch:cfg_util:savejob:unknown', 'Unknown file format for %s.', varargin{2});
            end
        end
    case 'setdef',
        % Set defaults for new jobs only
        cm = subsref(c0, [varargin{1}, varargin{2}]);
        cm = setval(cm, varargin{3}, true);
        c0 = subsasgn(c0, [varargin{1}, varargin{2}], cm);
    case 'setval',
        cjob        = varargin{1};
        mod_job_id  = varargin{2};
        item_mod_id = varargin{3};
        if cfg_util('ismod_job_id', cjob, mod_job_id) && ...
                cfg_util('isitem_mod_id', item_mod_id)
            cm = subsref(jobs(cjob).cj, [jobs(cjob).cjid2subs{mod_job_id}, item_mod_id]);
            cm = setval(cm, varargin{4}, false);
            jobs(cjob).cj = subsasgn(jobs(cjob).cj, [jobs(cjob).cjid2subs{mod_job_id}, item_mod_id], cm);
            varargout{1} = all_set_item(cm);
            % clear run configuration
            jobs(cjob).cjrun = [];
            jobs(cjob).cjid2subsrun = {};
        end
    case 'showdoc',
        if nargin == 2
            % get item from defaults tree
            cm = local_getcm(c0, varargin{1});
        elseif nargin >= 3
            % get item from job
            cm = local_getcmjob(jobs, varargin{:});
        end
        varargout{1} = showdoc(cm,'');
    case 'showdocwidth',
        if nargin == 3
            % get item from defaults tree
            cm = local_getcm(c0, varargin{2});
        elseif nargin >= 4
            % get item from job
            cm = local_getcmjob(jobs, varargin{2:end});
        end
        varargout{1} = cfg_justify(varargin{1}, showdoc(cm,''));
    case 'showjob',
        cjob = varargin{1};
        if cfg_util('isjob_id', cjob)
            if nargin > 2
                mod_job_ids = varargin{2};
                [unused str sts dep sout] = local_showjob(jobs(cjob).cj, ...
                                                          subsref(jobs(cjob).cjid2subs, ...
                                                                  substruct('{}', mod_job_ids)));
            else
                [id str sts dep sout] = local_showjob(jobs(cjob).cj, jobs(cjob).cjid2subs);
            end
            varargout{1} = id;
            varargout{2} = str;
            varargout{3} = sts;
            varargout{4} = dep;
            varargout{5} = sout;
        else
            varargout = {{}, {}, [], [], {}};
        end
    case 'tag2cfg_id',
        [mod_cfg_id item_mod_id] = local_tag2cfg_id(c0, varargin{1}, ...
                                                        true);
        if iscell(mod_cfg_id)
            % don't force mod_cfg_id to point to cfg_exbranch
            mod_cfg_id = local_tag2cfg_id(c0, varargin{1}, false);
        end
        varargout{1} = mod_cfg_id;
        varargout{2} = item_mod_id;
    case 'tag2mod_cfg_id',
        varargout{1} = local_tag2cfg_id(c0, varargin{1}, true);
    otherwise
        cfg_message('matlabbatch:usage', '%s: Unknown command ''%s''.', mfilename, cmd);
end
return;

%% Local functions
% These are the internal implementations of commands.
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [c0, jobs] = local_addapp(c0, jobs, cfg, varargin)
% Add configuration data to c0 and all jobs
% Input
% * cfg - Function name, function handle or cfg_item tree. If a
%         function is passed, this will be evaluated with no arguments and
%         must return a single configuration tree.
% * def - Optional. Function name, function handle or defaults struct/cell.
%         This function should return a job struct suitable to initialise
%         the defaults branches of the cfg tree.
%         If def is empty or does not exist, no defaults will be added.

if isempty(cfg)
    % Gracefully return if there is nothing to do
    return;
end
if subsasgn_check_funhandle(cfg)
    cvstate = cfg_get_defaults('cfg_item.checkval');
    cfg_get_defaults('cfg_item.checkval',true);
    c1 = feval(cfg);
    cfg_get_defaults('cfg_item.checkval',cvstate);
elseif isa(cfg, 'cfg_item')
    c1 = cfg;
end
if ~isa(c1, 'cfg_item')
    cfg_message('matlabbatch:cfg_util:addapp:inv',...
                'Invalid configuration');
    return;
end
dpind = cellfun(@(c0item)strcmp(c1.tag, c0item.tag), c0.values);
if any(dpind)
    dpind = find(dpind);
    cfg_message('matlabbatch:cfg_util:addapp:dup',...
                'Duplicate application tag in applications ''%s'' and ''%s''.', ...
                c1.name, c0.values{dpind(1)}.name);
    return;
end
if nargin > 3 && ~isempty(varargin{1})
    c1 = local_initdef(c1, varargin{1});
end
c0.values{end+1} = c1;
for k = 1:numel(jobs)
    jobs(k).cj.values{end+1} = c1;
    jobs(k).c0.values{end+1} = c1;
    % clear run configuration
    jobs(k).cjrun = [];
    jobs(k).cjid2subsrun = {};
end
cfg_message('matlabbatch:cfg_util:addapp:done', 'Added application ''%s''\n', c1.name);
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [job, id] = local_addtojob(job, c0subs)
% Add module subsref(c0, c0subs) to job.cj, append its subscript to
% job.cjid2subs and return the index into job.cjid2subs to the caller.
% The module will be added in a 'degenerated' branch of a cfg tree, where
% it is the only exbranch that can be reached on the 'val' path.
id = numel(job.cj.val)+1;
cjsubs = c0subs;
for k = 1:2:numel(cjsubs)
    % assume subs is [.val(ues){X}]+ and there are only choice/repeats
    % above exbranches
    % replace values{X} with val{1} in '.' references
    if strcmp(cjsubs(k).subs, 'values')
        cjsubs(k).subs = 'val';
        if k == 1
            % set id in cjsubs(2)
            cjsubs(k+1).subs = {id};
        else
            cjsubs(k+1).subs = {1};
        end
    end
    cm = subsref(job.c0, c0subs(1:(k+1)));
    if k == numel(cjsubs)-1
        % initialise dynamic defaults - not only in defaults tree, but
        % also in pre-configured .val items.
        cm = initialise(cm, '<DEFAULTS>', false);
    end
    % add path to module to cj
    job.cj = subsasgn(job.cj, cjsubs(1:(k+1)), cm);
end
% set id in module    
job.cj = subsasgn(job.cj, [cjsubs substruct('.', 'id')], cjsubs);
job.cjid2subs{id} = cjsubs;
% clear run configuration
job.cjrun = [];
job.cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function varargout = local_cd(pth)
% Try to work around some unexpected behaviour of MATLAB's cd command
if ~isempty(pth)
    if ischar(pth)
        wd = pwd;
        cd(pth);
    else
        cfg_message('matlabbatch:usage','CD: path must be a string.');
    end
else
    % Do not cd if pth is empty.
    wd = pwd;
end
if nargout > 0
    varargout{1} = wd;
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [job, n2oid] = local_compactjob(ojob)
% Remove placeholders from cj and recursively update dependencies to fit
% new ids. Warning: this will invalidate mod_job_ids!

job = ojob;
job.cj.val = {};
job.cjid2subs = {};

oid = cell(size(ojob.cjid2subs));
n2oid = NaN*ones(size(ojob.cjid2subs));
nid = 1;
for k = 1:numel(ojob.cjid2subs)
    if ~isempty(ojob.cjid2subs{k})
        cjsubs = ojob.cjid2subs{k};
        oid{nid} = ojob.cjid2subs{k};
        cjsubs(2).subs = {nid};
        job.cjid2subs{nid} = cjsubs;
        for l = 1:2:numel(cjsubs)
            % subs is [.val(ues){X}]+
            % add path to module to cj
            job.cj = subsasgn(job.cj, job.cjid2subs{nid}(1:(l+1)), ...
                          subsref(ojob.cj, ojob.cjid2subs{k}(1:(l+1))));
        end
        n2oid(nid) = k;
        nid = nid + 1;
    end
end
oid = oid(1:(nid-1));
n2oid = n2oid(1:(nid-1));
% update changed ids in job (where n2oid ~= 1:numel(n2oid))
cid = n2oid ~= 1:numel(n2oid);
if any(cid)
    job.cj = update_deps(job.cj, oid(cid), job.cjid2subs(cid));
end
% clear run configuration
job.cjrun = [];
job.cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function job = local_delfromjob(job, id)
% Remove module subsref(job.cj, job.cjid2subs{id}) from job.cj. All
% target and source dependencies between the module and other modules in
% job.cj are removed. Corresponding entries in job.cj and job.cjid2subs
% are set to {} in order to keep relationships within the tree consistent
% and in order to keep other ids valid. A rebuild of job.cj and an update
% of changed subsrefs would be possible (and needs to be done before
% e.g. saving the job). 
if isempty(job.cjid2subs) || isempty(job.cjid2subs{id}) || numel(job.cjid2subs) < id
    cfg_message('matlabbatch:cfg_util:invid', ...
            'Invalid id %d.', id);
    return;
end
cm = subsref(job.cj, job.cjid2subs{id});
if ~isempty(cm.tdeps)
    job.cj = del_in_source(cm.tdeps, job.cj);
end
if ~isempty(cm.sdeps)
    job.cj = del_in_target(cm.sdeps, job.cj);
end
% replace module with placeholder
cp = cfg_const;
cp.tag = 'deleted_item';
cp.val = {''};
cp.hidden = true;
% replace deleted module at top level, not at branch level
job.cj = subsasgn(job.cj, job.cjid2subs{id}(1:2), cp);
job.cjid2subs{id} = struct([]);
% clear run configuration
job.cjrun = [];
job.cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function local_gencode(c0, fname, tropts, preamble)
% Generate code, split at nodes matching stopspec (if stopspec is not
% empty). fname will be overwritten if tropts is empty (i.e. for single
% file output or subtrees). Note that some manual fixes may be required
% (function handles, variable/function names).
% If a preamble is passed as cellstr, it will be prepended to each
% generated file after the function... line. If no preamble is specified and
% a file <fname>_mlb_preamble.m exists in the folder where the
% configuration is being written, this file will be read and included
% literally.
if isempty(tropts)||isequal(tropts,cfg_tropts({{}},1,Inf,1,Inf,true)) || ...
        isequal(tropts, cfg_tropts({{}},1,Inf,1,Inf,false))
    tropts(1).clvl = 1;
    tropts(1).mlvl = Inf;
    tropts(1).cnt  = 1;
    [p funcname e] = fileparts(fname);
    [cstr tag] = gencode_item(c0, '', {}, [funcname '_'], tropts);
    funcname = [funcname '_' tag];
    fname = fullfile(p, [funcname '.m']);
    unpostfix = '';
    while exist(fname, 'file')
        cfg_message('matlabbatch:cfg_util:gencode:fileexist', ...
                ['While generating code for cfg_item: ''%s'', %s. File ' ...
                 '''%s'' already exists. Trying new filename - you will ' ...
                 'need to adjust generated code.'], ...
                c0.name, tag, fname);
        unpostfix = [unpostfix '1'];
        fname = fullfile(p, [funcname unpostfix '.m']);
    end
    fid = fopen(fname, 'wt');
    fprintf(fid, 'function %s = %s\n', tag, funcname);
    fprintf(fid, '%s\n', preamble{:});
    fprintf(fid, '%s\n', cstr{:});
    fclose(fid);
else
    % generate root level code
    [p funcname e] = fileparts(fname);
    [cstr tag] = gencode_item(c0, 'jobs', {}, [funcname '_'], tropts);
    fname = fullfile(p, [funcname '.m']);
    if nargin < 4 || isempty(preamble) || ~iscellstr(preamble)
        try
            fid = fopen(fullfile(p, [funcname '_mlb_preamble.m']),'r');
            ptmp = textscan(fid,'%s','Delimiter',sprintf('\n'));
            fclose(fid);
            preamble = ptmp{1};
        catch
            preamble = {};
        end
    end
    fid = fopen(fname, 'wt');
    fprintf(fid, 'function %s = %s\n', tag, funcname);
    fprintf(fid, '%s\n', preamble{:});
    fprintf(fid, '%s\n', cstr{:});
    fclose(fid);
    % generate subtree code - find nodes one level below stop spec
    tropts.mlvl = tropts.mlvl+1;
    [ids stop] = list(c0, tropts.stopspec, tropts);
    ids = ids(stop); % generate code for stop items only
    ctropts = cfg_tropts({{}},1,Inf,1,Inf,tropts.dflag);
    for k = 1:numel(ids)
        if ~isempty(ids{k}) % don't generate root level code again
            local_gencode(subsref(c0, ids{k}), fname, ctropts, preamble);
        end
    end
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [cj, cjid2subs] = local_getcjid2subs(cjin)
% Find ids of exbranches. 
% find ids
exspec = cfg_findspec({{'class', 'cfg_exbranch'}});
tropts = cfg_tropts({{'class', 'cfg_exbranch'}}, 1, Inf, 0, Inf, false);
cjid2subsin = list(cjin, exspec, tropts);
cjid2subs = cjid2subsin;
cj = cjin;
cancjid2subs = false(size(cjid2subsin));
for k = 1:numel(cjid2subs)
    % assume 1-D subscripts into val, and there should be at least one
    tmpsubs = [cjid2subs{k}(2:2:end).subs];
    cancjid2subs(k) = all(cellfun(@(cs)isequal(cs,1), tmpsubs(2:end)));
end
if all(cancjid2subs)
    idsubs = substruct('.', 'id');
    for k = 1:numel(cjid2subs)
        cj = subsasgn(cj, [cjid2subs{k} idsubs], ...
                      cjid2subs{k});
    end
else
    cj.val = cell(size(cjid2subs));
    idsubs = substruct('.', 'id');
    for k = 1:numel(cjid2subs)
        if cancjid2subs(k)
            % add path to module to cj
            cpath = subsref(cjin, cjid2subs{k}(1:2));
            cpath = subsasgn(cpath, [cjid2subs{k}(3:end) ...
                                idsubs], cjid2subs{k});
            cj = subsasgn(cj, cjid2subs{k}(1:2), cpath);
        else
            % canonicalise SPM5 batches to cj.val{X}.val{1}....val{1}
            % This would break dependencies, but in SPM5 batches there should not be any
            % assume subs is [.val{X}]+ and there are only choice/repeats
            % above exbranches
            for l = 2:2:numel(cjid2subs{k})
                if l == 2
                    cjid2subs{k}(l).subs = {k};
                else
                    cjid2subs{k}(l).subs = {1};
                end
                % add path to module to cj
                cpath = subsref(cjin, cjid2subsin{k}(1:l));
                % clear val field for nodes below exbranch
                if l < numel(cjid2subs{k})
                    cpath.val = {};
                else
                    cpath.id = cjid2subs{k};
                end
                cj = subsasgn(cj, cjid2subs{k}(1:l), cpath);
            end
        end
    end
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [cm, cfg_id] = local_getcm(c0, cfg_id)
if cfg_util('ismod_cfg_id', cfg_id)
    % This should better test something like 'iscfg_id'
    % do nothing
else
    [mod_cfg_id, item_mod_id] = cfg_util('tag2cfg_id', cfg_id);
    if isempty(mod_cfg_id)
        cfg_message('matlabbatch:cfg_util:invid', ...
              'Item with tag ''%s'' not found.', cfg_id);
    else
        cfg_id = [mod_cfg_id, item_mod_id];
        cm = subsref(c0, cfg_id);
    end
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function cm = local_getcmjob(jobs, job_id, mod_job_id, item_mod_id)
if nargin < 4
    item_mod_id = struct('subs',{},'type',{});
end
if cfg_util('isjob_id', job_id) && cfg_util('ismod_job_id', mod_job_id) ...
        && cfg_util('isitem_mod_id', item_mod_id)
    cm = subsref(jobs(job_id).cj, ...
                 [jobs(job_id).cjid2subs{mod_job_id} item_mod_id]);
else
    cfg_message('matlabbatch:cfg_util:invid', ...
          'Item not found.');
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function local_initapps
% add application data
if isdeployed
    cfg_mlbatch_appcfg_master;
else
    appcfgs = which('cfg_mlbatch_appcfg','-all');
    cwd = pwd;
    dirs = cell(size(appcfgs));
    for k = 1:numel(appcfgs)
        % cd into directory containing config file
        [p n e] = fileparts(appcfgs{k});
        local_cd(p);
        % try to work around MATLAB bug in symlink handling
        % only add application if this directory has not been visited yet
        dirs{k} = pwd;
        if ~any(strcmp(dirs{k}, dirs(1:k-1)))
            try
                [cfg def] = feval('cfg_mlbatch_appcfg');
                ests = true;
            catch
                ests = false;
                estr = cfg_disp_error(lasterror);
                cfg_message('matlabbatch:cfg_util:eval_appcfg', ...
                    'Failed to load %s', which('cfg_mlbatch_appcfg'));
                cfg_message('matlabbatch:cfg_util:eval_appcfg', '%s\n', estr{:});
            end
            if ests
                cfg_util('addapp', cfg, def);
            end
        end
    end
    local_cd(cwd);
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [c0, jobs, cjob] = local_initcfg
% initial config
c0   = cfg_mlbatch_root;
cjob = 1;
jobs(cjob).cj        = c0;
jobs(cjob).c0        = c0;
jobs(cjob).cjid2subs = {};
jobs(cjob).cjrun     = [];
jobs(cjob).cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function c1 = local_initdef(c1, varargin)
if nargin > 1
    defspec = varargin{1};
    if subsasgn_check_funhandle(defspec)
        opwd = pwd;
        if ischar(defspec)
            [p fn e] = fileparts(defspec);
            local_cd(p);
            defspec = fn;
        end
        def = feval(defspec);
        local_cd(opwd);
    elseif isa(defspec, 'cell') || isa(defspec, 'struct')
        def = defspec;
    else
        def = [];
    end
    if ~isempty(def)
        c1 = initialise(c1, def, true);
    end
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [cjob, mod_job_idlist] = local_initjob(cjob, job, jobdedup)
% Initialise a cell array of jobs
idsubs = substruct('.','id');
sdsubs = substruct('.','sdeps');
tdsubs = substruct('.','tdeps');
v1subs = substruct('.','val','{}',{1});
% Update application defaults - not only in .values, but also in
% pre-configured .val items. Use this as template for all jobs.
cjd = initialise(cjob.c0, '<DEFAULTS>', false);
% Prepare all jobs - initialise similar jobs only once.
[ujobdedup ui uj] = unique(jobdedup);
% Initialise jobs
ucj = cellfun(@(ucjob)initialise(cjd,ucjob,false), job(ui), 'UniformOutput', false);
% Canonicalise (this may break dependencies, see comment in
% local_getcjid2subs)
[ucj ucjid2subs] = cellfun(@local_getcjid2subs, ucj, 'UniformOutput', false);
% Harvest, keeping dependencies
[u1 u2 u3 u4 u5 ucj] = cellfun(@(cucj)harvest(cucj, cucj, false, false), ucj, 'UniformOutput', false);
% Deduplicate, concatenate
cjob.cj = ucj{uj(1)};
cjob.cjid2subs = ucjid2subs{uj(1)};
for n = 2:numel(uj)
    cjidoffset = numel(cjob.cjid2subs);
    cj1        = ucj{uj(n)};
    cjid2subs1 = ucjid2subs{uj(n)};
    cjid2subs2 = ucjid2subs{uj(n)};
    for k = 1:numel(cjid2subs2)
        % update id subscripts
        cjid2subs2{k}(2).subs{1} = cjid2subs2{k}(2).subs{1} + cjidoffset;
        cj1 = subsasgn(cj1, [cjid2subs1{k}, idsubs], ...
            cjid2subs2{k});
        % update src_exbranch in dependent cfg_items
        sdeps = subsref(cj1, [cjid2subs1{k}, sdsubs]);
        for l = 1:numel(sdeps)
            csdep = sdeps(l);
            tgt_exbranch = csdep.tgt_exbranch;
            tgt_input    = [csdep.tgt_input v1subs];
            % dependencies in dependent item
            ideps = subsref(cj1, [tgt_exbranch tgt_input]);
            isel = cellfun(@(iid)isequal(iid,cjid2subs1{k}),{ideps.src_exbranch});
            ideps(isel).src_exbranch = cjid2subs2{k};
            % save updated item & module
            cj1 = subsasgn(cj1, [tgt_exbranch tgt_input], ...
                ideps);
            % delete old tdeps - needs to be updated by harvest
            cj1 = subsasgn(cj1, [tgt_exbranch tdsubs], []);
        end
        % done with sdeps - clear
        cj1 = subsasgn(cj1, [cjid2subs1{k}, sdsubs], []);
    end
    % concatenate configs
    cjob.cjid2subs = [cjob.cjid2subs cjid2subs2];
    cjob.cj.val    = [cjob.cj.val cj1.val];
end
% harvest, update dependencies
[u1 u2 u3 u4 u5 cjob.cj] = harvest(cjob.cj, cjob.cj, false, false);
mod_job_idlist = num2cell(1:numel(cjob.cjid2subs));
% add field to keep run results from job
cjob.cjrun = [];
cjob.cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [job, id] = local_replmod(job, oid)
% Replicate module subsref(job.cj,job.cjid2subs{oid}) by adding it to the end of
% the job list. Update id in module and delete links to dependent modules,
% these dependencies are the ones of the original module, not of the
% replica.
id = numel(job.cj.val)+1;
% subsref of original module
ocjsubs = job.cjid2subs{oid};
% subsref of replica module
rcjsubs = ocjsubs;
rcjsubs(2).subs = {id};
for k = 1:2:numel(ocjsubs)
    % Add path to replica module, copying items from original path
    job.cj = subsasgn(job.cj, rcjsubs(1:(k+1)), subsref(job.cj, ocjsubs(1:(k+1))));
end
% set id in module, delete copied sdeps and tdeps
job.cj = subsasgn(job.cj, [rcjsubs substruct('.', 'id')], rcjsubs);
job.cj = subsasgn(job.cj, [rcjsubs substruct('.', 'sdeps')], []);
job.cj = subsasgn(job.cj, [rcjsubs substruct('.', 'tdeps')], []);
% re-harvest to update tdeps and outputs
[u1 u2 u3 u4 u5 job.cj] = harvest(subsref(job.cj, rcjsubs), job.cj, false, false);
job.cjid2subs{id} = rcjsubs;
% clear run configuration
job.cjrun = [];
job.cjid2subsrun = {};
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [job err] = local_runcj(job, cjob, pflag, cflag)
% Matlab uses a copy-on-write policy with very high granularity - if
% modified, only parts of a struct or cell array are copied.
% However, forward resolution may lead to high memory consumption if
% variables are passed, but avoids extra housekeeping for outputs and
% resolved dependencies.
% Here, backward resolution is used. This may be time consuming for large
% jobs with many dependencies, because dependencies of any cfg_item are
% resolved only if all of them are resolvable (i.e. there can't be a mix of
% values and dependencies in a .val field).
% If pflag is true, then modules will be executed in parallel if they are
% independent. Setting pflag to false forces serial execution of modules
% even if they seem to be independent.
% If cflag is true, and a job with pre-set module outputs .jout is passed
% in job.cjrun, the corresponding modules will not be run again.

cfg_message('matlabbatch:run:jobstart', ...
            ['\n\n------------------------------------------------------------------------\n',...
             'Running job #%d\n', ...
             '------------------------------------------------------------------------'], cjob);
if cflag && ~isempty(job.cjrun)
    [u1 mlbch] = harvest(job.cjrun, job.cjrun, false, true);
else
    job1 = local_compactjob(job);
    job.cjid2subsrun = job1.cjid2subs;
    [u1 mlbch u3 u4 u5 job.cjrun] = harvest(job1.cj, job1.cj, false, true);
end
% copy cjid2subs, it will be modified for each module that is run
cjid2subs = job.cjid2subsrun;
cjid2subsfailed = {};
cjid2subsskipped = {};
tdsubs = substruct('.','tdeps');
chsubs = substruct('.','chk');
while ~isempty(cjid2subs)
    % find mlbch that can run
    cand = false(size(cjid2subs));
    if pflag
        % Check dependencies of all remaining mlbch
        maxcand = numel(cjid2subs);
    else
        % Check dependencies of first remaining job only
        maxcand = min(1, numel(cjid2subs));
    end
    for k = 1:maxcand
        cand(k) = isempty(subsref(job.cjrun, [cjid2subs{k} tdsubs])) ...
                  && subsref(job.cjrun, [cjid2subs{k} chsubs]) ...
                  && all_set(subsref(job.cjrun, cjid2subs{k}));
    end
    if ~any(cand)
        cfg_message('matlabbatch:run:nomods', ...
                'No executable modules, but still unresolved dependencies or incomplete module inputs.');
        cjid2subsskipped = cjid2subs;
        break;
    end
    % split job list
    cjid2subsrun = cjid2subs(cand);
    cjid2subs = cjid2subs(~cand);
    % collect sdeps of really running modules
    csdeps = cell(size(cjid2subsrun));
    % run modules that have all dependencies resolved
    for k = 1:numel(cjid2subsrun)
        cm = subsref(job.cjrun, cjid2subsrun{k});
        if isa(cm.jout,'cfg_inv_out')
            % no cached outputs (module did not run or it does not return
            % outputs) - run job
            cfg_message('matlabbatch:run:modstart', 'Running ''%s''', cm.name);
            try
                cm = cfg_run_cm(cm, subsref(mlbch, cfg2jobsubs(job.cjrun, cjid2subsrun{k})));
                csdeps{k} = cm.sdeps;
                cfg_message('matlabbatch:run:moddone', 'Done    ''%s''', cm.name);
            catch
                cjid2subsfailed = [cjid2subsfailed cjid2subsrun(k)];
                le = lasterror;
                % try to filter out stack trace into matlabbatch
                try
                    runind = find(strcmp('cfg_run_cm', {le.stack.name}));
                    le.stack = le.stack(1:runind-1);
                end
                str = cfg_disp_error(le);
                cfg_message('matlabbatch:run:modfailed', 'Failed  ''%s''', cm.name);
                cfg_message('matlabbatch:run:modfailed', '%s\n', str{:});
            end
            % save results (if any) into job tree
            job.cjrun = subsasgn(job.cjrun, cjid2subsrun{k}, cm);
        else
            % Use cached outputs
            cfg_message('matlabbatch:run:cached', 'Using cached outputs for ''%s''', cm.name);
        end
    end
    % update dependencies, re-harvest mlbch
    tmp = [csdeps{:}];
    if ~isempty(tmp)
        ctgt_exbranch = {tmp.tgt_exbranch};
        % assume job.cjrun.val{k}.val{1}... structure
        ctgt_exbranch_id = zeros(size(ctgt_exbranch));
        for k = 1:numel(ctgt_exbranch)
            ctgt_exbranch_id(k) = ctgt_exbranch{k}(2).subs{1};
        end
        % harvest only targets and only once
        [un ind] = unique(ctgt_exbranch_id);
        for k = 1:numel(ind)
            cm = subsref(job.cjrun, ctgt_exbranch{ind(k)});
            [u1 cmlbch u3 u4 u5 job.cjrun] = harvest(cm, job.cjrun, false, ...
                                                     true);
            mlbch = subsasgn(mlbch, cfg2jobsubs(job.cjrun, ctgt_exbranch{ind(k)}), ...
                             cmlbch);
        end
    end
end
if isempty(cjid2subsfailed) && isempty(cjid2subsskipped)
    cfg_message('matlabbatch:run:jobdone', 'Done\n');
    err = [];
else
    str = cell(numel(cjid2subsfailed)+numel(cjid2subsskipped)+1,1);
    str{1} = 'The following modules did not run:';
    for k = 1:numel(cjid2subsfailed)
        str{k+1} = sprintf('Failed: %s', subsref(job.cjrun, [cjid2subsfailed{k} substruct('.','name')]));
    end
    for k = 1:numel(cjid2subsskipped)
        str{numel(cjid2subsfailed)+k+1} = sprintf('Skipped: %s', ...
                                                  subsref(job.cjrun, [cjid2subsskipped{k} substruct('.','name')]));
    end
% Commented out for SPM public release
%    str{end+1} = sprintf(['If the problem can be fixed without modifying ' ...
%                        'the job, the computation can be resumed by ' ...
%                        'running\n   cfg_util(''cont'',%d)\nfrom the ' ...
%                        'MATLAB command line.'],cjob);
    cfg_message('matlabbatch:run:jobfailed', '%s\n', str{:});
    err.identifier = 'matlabbatch:run:jobfailederr';
    err.message    = sprintf(['Job execution failed. The full log of this run can ' ...
                        'be found in MATLAB command window, starting with ' ...
                        'the lines (look for the line showing the exact ' ...
                        '#job as displayed in this error message)\n' ...
                        '------------------ \nRunning job #%d' ...
                        '\n------------------\n'], cjob);
    err.stack      = struct('file','','name','MATLABbatch system','line',0);
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [id, str, sts, dep, sout] = local_showjob(cj, cjid2subs)
% Return name, all_set status and id of internal job representation
id  = cell(size(cjid2subs));
str = cell(size(cjid2subs));
sts = false(size(cjid2subs));
dep = false(size(cjid2subs));
sout = cell(size(cjid2subs));
cmod = 1; % current module count
for k = 1:numel(cjid2subs)
    if ~isempty(cjid2subs{k})
        cm = subsref(cj, cjid2subs{k});
        id{cmod}  = k;
        str{cmod} = cm.name;
        sts(cmod) = all_set(cm);
        dep(cmod) = ~isempty(cm.tdeps);
        sout{cmod} = cm.sout;
        cmod = cmod + 1;
    end
end
id   = id(1:(cmod-1));
str  = str(1:(cmod-1));
sts  = sts(1:(cmod-1));
dep  = dep(1:(cmod-1));
sout = sout(1:(cmod-1));
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function [mod_cfg_id, item_mod_id] = local_tag2cfg_id(c0, tagstr, splitspec)

tags = textscan(tagstr, '%s', 'delimiter', '.');
taglist = tags{1};
if ~strcmp(taglist{1}, c0.tag)
    % assume tag list starting at application level
    taglist = [c0.tag taglist(:)'];
end
if splitspec
    % split ids at cfg_exbranch level
    finalspec = cfg_findspec({{'class','cfg_exbranch'}});
else
    finalspec = {};
end
tropts=cfg_tropts({{'class','cfg_exbranch'}},0, inf, 0, inf, true);
[mod_cfg_id stop rtaglist] = tag2cfgsubs(c0, taglist, finalspec, tropts);
if iscell(mod_cfg_id)
    item_mod_id = {};
    return;
end

if isempty(rtaglist)
    item_mod_id = struct('type',{}, 'subs',{});
else
    % re-add tag of stopped node
    taglist = [gettag(subsref(c0, mod_cfg_id)) rtaglist(:)'];
    tropts.stopspec = {};
    [item_mod_id stop rtaglist] = tag2cfgsubs(subsref(c0, mod_cfg_id), ...
                                              taglist, {}, tropts);
end
