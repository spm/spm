function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)

% function [tag, val, typ, dep, chk, cj] = harvest(item, cj, dflag, rflag)
% Harvest a cfg_branch object.
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
% $Id: harvest.m 2512 2008-12-01 13:21:29Z volkmar $

rev = '$Rev: 2512 $'; %#ok

typ = class(item);
tag = gettag(item);
val = struct([]);
dep = []; % Placeholder for dependencies. Will be classified during
          % first call to dep_add
chk = ~dflag && rflag;

tname = treepart(item, dflag);
ntgt_input = substruct('.', tname, '{}', {});
citems = subsref(item, ntgt_input(1));

% add references into harvested struct/cell
njtsubs.type = '.';
njtsubs.subs = '';

for k = 1:numel(citems)
    [ctag cval unused cdep cchk cj] = harvest(citems{k}, cj, dflag, rflag);
    val(1).(ctag) = cval;
    if ~dflag && ~isempty(cdep)
        % augment cdep tsubs references
        njtsubs.subs  = ctag;
        ntgt_input(2).subs  = {k};
        dep = dep_add(cdep, dep, ntgt_input, njtsubs);
    end;
    chk = chk && cchk;
end;
if chk
    chk = docheck(item, val);
end;
