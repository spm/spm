function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)

% function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)
% Generic harvest function, suitable for all const/entry items.
% The configuration tree cj is passed unmodified. If rflag is true and a
% dependency can be resolved, the resolved value will be returned,
% otherwise the cfg_dep object will be returned in val and dep.
% Input arguments:
% item  - item to be harvested
% cj    - configuration tree (passed unmodified)
% dflag - if true, harvest defaults tree, otherwise filled tree
% rflag - if true, resolve dependencies in leaf nodes
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
% $Id: harvest.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

typ = class(item);
tag = item.tag;
val = '<UNDEFINED>';
dep = []; % placeholder for dependencies
chk = ~dflag && rflag;

if ~isempty(item.val)
    if isa(item.val{1},'cfg_dep')
        if dflag % do not harvest references if defaults are requested
            val = [];
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
if chk && ~isempty(item.check)
    % Call content dependent check function, if it is present and
    % dependencies are to be resolved
    try
        cstr = feval(item.check, val);
        chk  = isempty(cstr);
    catch
        cstr = sprintf('Check function ''%s'' failed.', ...
                       func2str(item.check));
        chk  = false;
    end;
    if ~chk
        warning('matlabbatch:harvest:check', ...
                'Contents of ''%s'' does not meet check criteria: ''%s''', ...
                tag, cstr);
    end;
end;
