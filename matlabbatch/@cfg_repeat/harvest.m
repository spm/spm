function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)

% function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)
% Harvest a cfg_repeat object.
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
% $Id: harvest.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

typ = class(item);
tag = gettag(item);
dep = []; % Placeholder for dependencies. Will be classified during
          % first call to dep_add.
chk = ~dflag && rflag;

tname = treepart(item, dflag);
ntgt_input = substruct('.', tname, '{}', {});
citems = subsref(item, ntgt_input(1));
if numel(item.values)==1 && isa(item.values{1},'cfg_branch') && ~item.forcestruct,
    if numel(citems) == 0
        % return empty struct
        cargs = cell(2*numel(item.values{1}.val),1);
        for i=1:numel(item.values{1}.val),
            cargs{2*i-1} = gettag(item.values{1}.val{i});
            cargs{2*i}   = {};
        end;
        val = struct(cargs{:});
    end;
    if ~dflag
        njtsubs(1).type = '()';
    end;
else
    if ~dflag
        val = cell(size(citems));
        njtsubs(1).type = '{}';
    else
        val = struct([]);
    end;
end;
for i=1:numel(citems),
    [ctag, cval, unused, cdep, cchk, cj] = harvest(citems{i}, cj, dflag, rflag);
    if numel(item.values)==1 && isa(item.values{1},'cfg_branch') && ~item.forcestruct,
        % FIXME: don't know how to best preinit this without raising
        % warnings about incompatible assignments ...
        val(i) = cval;
    else
        if numel(item.values)>1 || item.forcestruct,
            if iscell(cval)
                cval = struct(ctag,{cval});
            else
                cval = struct(ctag,cval);
            end;
            if ~dflag
                njtsubs(2).type = '.';
                njtsubs(2).subs  = ctag;
            end;
        end;
        if dflag
            % return a struct containing defaults of all child nodes
            % instead of a cell array. This makes defaults easier to
            % read.
            if numel(item.values)>1 || item.forcestruct,
                val(1).(ctag) = cval.(ctag);
            else
                val = {cval};
            end;
        else
            val{i} = cval;
        end;
    end;
    if ~dflag && ~isempty(cdep)
        % augment cdep tsubs references
        ntgt_input(2).subs = {i};
        njtsubs(1).subs = {i};
        dep = dep_add(cdep, dep, ntgt_input, njtsubs);
    end;
    chk = chk && cchk;
end;
if chk 
    chk = docheck(item, val);
end;
