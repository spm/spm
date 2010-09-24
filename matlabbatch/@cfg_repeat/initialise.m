function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Initialise a configuration tree with values. If val is a job
% struct/cell, only the parts of the configuration that are present in
% this job will be initialised. If dflag is true, then matching items
% from item.values will be initialised. If dflag is false, matching items
% from item.values will be added to item.val and initialised after
% copying.
% If val has the special value '<DEFAULTS>', the entire configuration
% will be updated with values from .def fields. If a .def field is
% present in a cfg_leaf item, the current default value will be inserted,
% possibly replacing a previously entered (default) value. If dflag is
% true, defaults will only be set in item.values. If dflag is false,
% defaults will be set for both item.val and item.values.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 4073 2010-09-24 12:07:57Z volkmar $

rev = '$Rev: 4073 $'; %#ok

if strcmp(val,'<DEFAULTS>')
    item = initialise_def(item, val, dflag);
else
    item = initialise_job(item, val, dflag);
end;

function item = initialise_def(item, val, dflag)
if ~dflag
    % initialise defaults both in current job and in defaults
    citem = subsref(item, substruct('.','val'));
    for k = 1:numel(citem)
        citem{k} = initialise(citem{k}, val, dflag);
    end;
    item = subsasgn(item, substruct('.','val'), citem);
end;
for k = 1:numel(item.values)
    item.values{k} = initialise(item.values{k}, val, dflag);
end;

function item = initialise_job(item, val, dflag)
if numel(item.values)==1 && isa(item.values{1},'cfg_branch') ...
        && ~item.forcestruct,
    if isstruct(val)
        if dflag
            item.values{1} = initialise(item.values{1}, val, dflag);
        else
            citem = cell(1,numel(val));
            for k = 1:numel(val)
                citem{k} = initialise(item.values{1}, val(k), dflag);
            end;
            item.cfg_item.val = citem;
        end;
    else
        cfg_message('matlabbatch:initialise', ...
                    'Can not initialise %s value(s): job is not a struct.', ...
                    gettag(item));
        return;
    end;
else
    if dflag
        if numel(item.values) > 1 || item.forcestruct
            % val should be either a cell array containing structs with a
            % single field (a harvested job), or a struct with multiple
            % fields (a harvested defaults tree). In the latter case,
            % convert val to a cell array before proceeding.
            if isstruct(val)
                vtag = fieldnames(val);
                val1 = cell(size(vtag));
                for k = 1:numel(vtag)
                    val1{k} = struct(vtag{k}, {val.(vtag{k})});
                end;
                val = val1;
            elseif iscell(val) && all(cellfun(@isstruct,val))
                vtag = cell(size(val));
                for k = 1:numel(val)
                    vtag(k) = fieldnames(val{k});
                end;
            else
                cfg_message('matlabbatch:initialise', ...
                            'Can not initialise %s value(s): job is not a cell array of struct items.', ...
                            gettag(item));
                return;
            end;
            for k = 1:numel(item.values)
                % use first match for defaults initialisation
                sel = find(strcmp(gettag(item.values{k}), vtag));
                if ~isempty(sel)
                    item.values{k} = initialise(item.values{k}, ...
                                                val{sel(1)}.(vtag{sel(1)}), ...
                                                dflag);
                end;
            end;
        else
            item.values{1} = initialise(item.values{1}, val{1}, dflag);
        end;
    else
        citem = cell(1,numel(val));
        if numel(item.values) > 1 || item.forcestruct
            for l = 1:numel(val)
                % val{l} should be a struct with a single field
                vtag = fieldnames(val{l});
                for k = 1:numel(item.values)
                    if strcmp(gettag(item.values{k}), vtag{1})
                        citem{l} = initialise(item.values{k}, ...
                            val{l}.(vtag{1}), ...
                            dflag);
                    end;
                end;
            end;
        else
            for l = 1:numel(val)
                citem{l} = initialise(item.values{1}, ...
                    val{l}, dflag);
            end;
        end;
        item.cfg_item.val = citem;
    end;
end;

