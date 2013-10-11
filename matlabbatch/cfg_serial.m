function cfg_serial(guifcn, job, varargin)

% This function is deprecated.
% The functionality should replaced by the following sequence of calls:
%
% Instead of
% cfg_serial(guifcn, job, varargin)
% use
% cjob = cfg_util('initjob', job);
% sts  = cfg_util('filljobui', cjob, guifcn, varargin);
% if sts
%      cfg_util('run', cjob);
% end;
% cfg_util('deljob', cjob);
%
% Instead of
% cfg_serial(guifcn, tagstr, varargin)
% use
% cjob = cfg_util('initjob');
% mod_cfg_id = cfg_util('tag2cfg_id', tagstr);
% cfg_util('addtojob', cjob, mod_cfg_id);
% sts  = cfg_util('filljobui', cjob, guifcn, varargin);
% if sts
%      cfg_util('run', cjob);
% end;
% cfg_util('deljob', cjob);
%
% Instead of
% cfg_serial(guifcn, mod_cfg_id, varargin)
% use
% cjob = cfg_util('initjob');
% cfg_util('addtojob', cjob, mod_cfg_id);
% sts  = cfg_util('filljobui', cjob, guifcn, varargin);
% if sts
%      cfg_util('run', cjob);
% end;
% cfg_util('deljob', cjob);
%
% If no guifcn is specified, use cfg_util('filljob',... instead.
%
% GuiFcn semantics
% [val sts] = guifcn(item)
% val should be suitable to set item.val{1} using setval(item, val,
% false) for all cfg_leaf items. For cfg_repeat/cfg_choice items, val
% should be a cell array of indices into item.values. For each element of
% val, setval(item, [val{k} Inf], false)
% will be called and thus item.values{k} will be appended to item.val.
% sts should be set to true, if guifcn returns with success (i.e. a
% valid value is returned or input should continue for the next item,
% regardless of value validity).

% Old help
% function cfg_serial(guifcn, job|tagstr|mod_cfg_id, varargin)
% A matlabbatch user interface which completes a matlabbatch job with
% incomplete inputs one at a time and runs the job.
% This interface may be called with or without user interaction. If
% guifcn is a function handle, then this function will be called to enter
% unspecified inputs. Otherwise, the first argument will be ignored and
% inputs will be assigned from the argument list only.
% During job completion and execution, this interface may interfere with
% cfg_ui. However, after the job is executed, it will be removed from the
% job list.
%
% cfg_serial(guifcn, job|tagstr|mod_cfg_id[, input1, input2, ...inputN])
% Ask for missing inputs in a job. Job should be a matlabbatch job as
% returned by cfg_util('harvest'). Modifications to the job structure are
% limited:
% - no new modules can be added
% - no dependencies can be created
% - unset cfg_choice items can be selected, but not altered
% - unset cfg_repeat items can be completed, but no repeats can be added
%   to or removed from already set items
% Data to complete a job can specified as additional arguments. This
% allows to script filling of a predefined job with missing inputs.
% For cfg_entry and cfg_file items, this can be any data that meets the 
% consistency (type, filter, size) constraints of the item.
% For cfg_menu and cfg_choice items, it can be the number of a menu
% option, counting from 1.
% For cfg_repeat items, it is a cell list of numbers of menu options. The
% corresponding option will be added at the end of the repeat list.
% Any input data that can not be assigned to the current input item will
% be discarded.
%
% cfg_serial(guifcn, tagstr)
% Start input at the node addressed by tagstr in the configuration
% tree. If tagstr points to an cfg_exbranch, then only this cfg_exbranch
% will be filled and run. If tagstr points to a node above cfg_exbranch
% level, then multiple cfg_exbranches below this node may be added and
% filled.
%
% guifcn Interface 
% The guifcn will be called to enter inputs to cfg_choice, cfg_repeat,
% cfg_entry, cfg_files, cfg_menu items. The general call syntax is
% val = feval(guifcn, class, name, item_specific_data);
% Input arguments
%   class - string: one of cfg_choice, cfg_repeat, 
%                        cfg_entry, cfg_files, cfg_menu
%   name  - string: Item display name
%   Item specific data
%   * cfg_choice, cfg_menu
%     labels - cell array of strings, menu labels
%     values - cell array of values corresponding to labels
%   * cfg_repeat
%     labels and values as above
%     num    - 2-vector of min/max number of required elements in val
%   * cfg_entry
%     strtype - input type
%     extras  - extra data to evaluate input
%     num     - required size of data
%   * cfg_files
%     filter, ufilter - file filter specifications (see cfg_select)
%     dir - initial directory
%     num - 2-vector of min/max number of required files
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_serial.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

cfg_message('matlabbatch:deprecated:cfg_serial', '''cfg_serial'' is deprecated. Please use cfg_util(''filljob[ui]'',...) to fill a job in serial mode.');
if ischar(job)
    % Assume dot delimited sequence of tags
    mod_cfg_id = cfg_util('tag2mod_cfg_id', job);
    if cfg_util('ismod_cfg_id', mod_cfg_id)
        % tag string points to somewhere in an cfg_exbranch
        % initialise new job
        cjob = cfg_util('initjob');
        % add mod_cfg_id
        cfg_util('addtojob', cjob, mod_cfg_id);
    else
        % tag string points to somewhere above cfg_branch
        cjob = local_addtojob(job);
    end;
elseif cfg_util('ismod_cfg_id', job)
        % initialise new job
        cjob = cfg_util('initjob');
        % add mod_cfg_id
        cfg_util('addtojob', cjob, job);    
else
    % assume job to be a saved job structure
    cjob = cfg_util('initjob', job);
end;
% varargin{:} is a list of input items
in = varargin;
% get job information
[mod_job_idlist, str, sts, dep, sout] = cfg_util('showjob', cjob);
for cm = 1:numel(mod_job_idlist)
    % loop over modules, enter missing inputs
    if ~sts(cm)
        in = local_fillmod(guifcn, cjob, mod_job_idlist{cm}, in);
    end;
end;
cfg_util('run',cjob);
cfg_util('deljob',cjob);

%% local functions

function cjob = local_addtojob(job)
% traverse tree down to cfg_exbranch level, add selected modules to job
cfg_message('matlabbatch:cfg_serial:notimplemented', ...
      'Menu traversal not yet implemented.');

function inputs = local_fillmod(guifcn, cjob, cm, inputs)
[item_mod_idlist, stop, contents] = ...
    cfg_util('listmod', cjob, cm, [], cfg_findspec({{'hidden',false}}), ...
             cfg_tropts({{'hidden', true}},1,Inf,1,Inf,false), ...
             {'class', 'all_set_item'});
for ci = 1:numel(item_mod_idlist)
    % loop over items, enter missing inputs
    if ~contents{2}{ci}
        if ~isempty(inputs)
            sts = local_setval(cjob, cm, item_mod_idlist, contents, ci, inputs{1});
            % discard input{1}, even if setval failed
            inputs = inputs(2:end);
        else
            sts = false;
        end;
        if ~sts && ~isa(guifcn, 'function_handle')
            % no input given, or input did not match required criteria
            cfg_message('matlabbatch:cfg_serial:notimplemented', ...
                  'User prompted input not yet implemented.');
        end;
        while ~sts
            % call guifcn until a valid input is returned
            val = local_call_guifcn(guifcn, cjob, cm, item_mod_idlist{ci}, ...
                                    contents{1}{ci});
            sts = local_setval(cjob, cm, item_mod_idlist, contents, ci, val);
        end;
        if strcmp(contents{1}{ci}, 'cfg_choice')||...
                strcmp(contents{1}{ci}, 'cfg_repeat')
            % restart filling current module, break out of for loop
            % afterwards
            inputs = local_fillmod(guifcn, cjob, cm, inputs);
            return;
        end;
    end;
end;

function sts = local_setval(cjob, cm, item_mod_idlist, contents, ci, val)
if strcmp(contents{1}{ci}, 'cfg_repeat')
    % assume val to be a cell array of indices into
    % .values
    % note that sts may return true even if some of the
    % indices failed. sts only indicates that the cfg_repeat
    % all_set_item status is met (i.e. the min/max number of
    % repeated items are present).
    sts = false;
    for cv = 1:numel(val)
        % do not use fast track || here, otherwise
        % cfg_util('setval') must be called in any case to
        % append val{cv} to cfg_repeat list.
        sts = sts | cfg_util('setval', cjob, cm, item_mod_idlist{ci}, ...
                             [val{cv} Inf]);
    end;
else
    % try to set val
    sts = cfg_util('setval', cjob, cm, item_mod_idlist{ci}, ...
                   val);
end;

function val = local_call_guifcn(guifcn, cjob, cm, citem, cmclass)
% fieldnames depend on class of item
switch cmclass
    case {'cfg_choice', 'cfg_repeat'},
        fnames = {'name', 'values', 'num'};
    case 'cfg_entry',
        fnames = {'name', 'strtype', 'num', 'extras'};
    case 'cfg_files',
        fnames = {'name', 'num', 'filter', 'dir', 'ufilter'};
    case 'cfg_menu',
        fnames = {'name', 'labels', 'values'};
end;

% only search current module/item
fspec  = cfg_findspec({{'hidden',false}});
tropts = cfg_tropts({{'hidden', true}},1,1,1,1,false);
if isempty(citem)
    % we are not in a module
    [u1, u2, contents] = cfg_util('listcfgall', cm, fspec, tropts, fnames);
else
    [u1, u2, contents] = cfg_util('listmod', cjob, cm, citem, fspec, tropts, fnames);
end;

cmname = contents{1}{1};
switch cmclass
    case {'cfg_choice', 'cfg_repeat'}
        % construct labels and values from 'values' field
        labels = cell(size(contents{2}{1}));
        for k = 1:numel(contents{2}{1})
            labels{k} = contents{2}{1}{k}.name;
        end;
        values = num2cell(1:numel(contents{2}{1}));
        args = {labels, values};
        if strcmp(cmclass, 'cfg_repeat')
            args{3} = contents{3}{1};
        end;
    case {'cfg_entry', 'cfg_files', 'cfg_menu'}
        args = cell(1, numel(contents)-1);
        for k = 2:numel(contents)
            args{k-1} = contents{k}{1};
        end;
end;
val = feval(guifcn, cmclass, cmname, args{:});
