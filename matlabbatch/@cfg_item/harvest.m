function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)

% function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)
% This is the generic harvest function, suitable for all cfg_leaf items.
% The configuration tree cj is passed unmodified. If rflag is true and a
% dependency can be resolved, the resolved value will be returned,
% otherwise the cfg_dep object will be returned in val and dep.
% If .val is empty and .def is set, the default value for this item will be
% returned instead.
% Input arguments:
% item  - item to be harvested
% cj    - configuration tree (passed unmodified)
% dflag - if true, harvest defaults tree, otherwise filled tree
% rflag - if true, resolve dependencies in cfg_leaf nodes
% Output arguments:
% tag - tag of harvested item
% val - harvested value
% typ - class of harvested item (currently unused)
% dep - list of unresolved dependencies
% chk - meaningful if ~dflag and all dependencies are resolved. Then it
%       returns success status of this items .check function and its
%       childrens check functions. A job is ready to run if all
%       dependencies are resolved and chk status is true.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: harvest.m 1517 2008-04-29 15:46:08Z volkmar $

rev = '$Rev: 1517 $';

typ = class(item);
tag = item.tag;
dep = cfg_dep;    % placeholder for dependencies
dep = dep(false); % make dep an empty dependency array
chk = ~dflag && rflag;

if isempty(item.val)
    if isempty(item.def)
        val = '<UNDEFINED>';
    else
        val = getdef(item);
    end;
else
    if isa(item.val{1},'cfg_dep')
        if dflag % do not harvest references if defaults are requested
            if isempty(item.def)
                val = '<UNDEFINED>';
            else
                val = getdef(item);
            end;
        else
            if rflag
                [val sts] = resolve_deps(item, cj);
            else
                sts = false;
            end;
            if ~sts % deps not resolved
                dep = item.val{1}; % dep is a array of cfg_dep objects
                for k = 1:numel(dep) % we may have multiple dependencies
                    dep(k).tname = item.name; % set target name
                end;
                val = dep; % return deps also in val for saving
            end;
        end;
    else
        val = item.val{1};
    end;
end;
chk = chk && isempty(dep);
if chk
    chk = docheck(item, val);
end;

function val = getdef(item)
try
    val = feval(item.def{:});
    if ~strcmp(val,'<UNDEFINED>')
        [sts val] = subsasgn_check(item, substruct('.','val'),{val});
        if sts
            % de-reference after subsasgn_check
            val = val{1};
        else
            val = '<UNDEFINED>';
        end;
    end;
catch
    val = '<UNDEFINED>';
    warning('matlabbatch:cfg_item:harvest:nodef', '%s: No matching defaults value found.', subsasgn_checkstr(item,substruct('.','val')));
end;
