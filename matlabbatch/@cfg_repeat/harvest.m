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
% $Id: harvest.m 1195 2008-03-07 21:51:49Z volkmar $

rev = '$Rev: 1195 $';

typ = class(item);
tag = gettag(item);
dep = [];
chk = ~dflag && rflag;

tname = treepart(item, dflag);
ntgt_input = substruct('.', tname, '{}', {});
citems = subsref(item, ntgt_input(1));
if numel(item.values)==1 && isa(item.values{1},'cfg_branch') && ~item.forcestruct,
    if numel(citems) == 0
        % initialise to empty struct
        cargs = {};
        for i=1:numel(item.values{1}.val),
            cargs = {cargs{:},gettag(item.values{1}.val{i}),{}};
        end;
        val = struct(cargs{:});
    end;
    if ~dflag
        njtsubs(1).type = '()';
    end;
else
    val = {};
    if ~dflag
        njtsubs(1).type = '{}';
    end;
end;
for i=1:numel(citems),
    [ctag cval unused cdep cchk cj] = harvest(citems{i}, cj, dflag, rflag);
    if numel(item.values)==1 && isa(item.values{1},'cfg_branch') && ~item.forcestruct,
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
                val.(ctag) = cval.(ctag);
            else
                val = {cval};
            end;
        else
            val = {val{:}, cval};
        end;
    end;
    if ~dflag
        ntgt_input(2).subs = {i};
        njtsubs(1).subs = {i};
        % augment cdep tsubs references
        dep = cfg_dep_add(dep, cdep, ntgt_input, njtsubs);
    end;
    chk = chk && cchk;
end;
if chk 
    chk = docheck(item, val);
end;
