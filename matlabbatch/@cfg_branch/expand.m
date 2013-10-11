function [item, sts] = expand(item, eflag, tropts)

% function [item, sts] = expand(item, eflag, tropts)
% Set/query expanded flag of item depending on eflag:
% -1 - do not force eflag to any state, only child state will be inherited
%  0 - collapse
%  1 - expand val unconditionally
%  2 - expand metadata unconditionally
%  3 - expand val, if it is not set
% Return status is (expanded > 0), i.e. if expanded, then no additional
% info about expansion level or expansion reason is returned and parent
% nodes are set to expanded = 1.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop traversal
% dflag    - traverse val or values tree
% clvl     - current level in tree
% mlvl     - maximum level to traverse - range 1 (top level only) to
%            Inf (all levels)
% cnt (not set here)
% mcnt (not evaluated here)
% Traversal options are used here to control which items should be forced
% to expand/unexpand. Traversal continues to child items, even if level or
% stopspec criteria are met, but with an eflag of -1 (i.e. only 'expanded'
% status is queried, but not changed).
%
% This function is identical for all cfg_intree classes.
%
% This code is part of a batch job configuration system for MATLAB. See
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: expand.m 5678 2013-10-11 14:58:04Z volkmar $

rev = '$Rev: 5678 $'; %#ok

% Set expanded based on eflag in item
if eflag >= 0 && eflag <= 2
    item.cfg_item.expanded = eflag;
end;
if eflag == 3
    if ~all_set_item(item)
        item.cfg_item.expanded = eflag;
    else
        item.cfg_item.expanded = 0;
    end;
end;
sts = item.cfg_item.expanded > 0;

% Check whether to proceed with forcing eflags
if eflag >= 0 && (tropts.clvl >= tropts.mlvl || (~isempty(tropts.stopspec) && match(item, tropts.stopspec)))
    eflag = -1;
end;

tname = treepart(item, tropts.dflag);
ntgt_input = substruct('.', tname);
citems = subsref(item, ntgt_input);

tropts.clvl = tropts.clvl + 1;

for k = 1:numel(citems)
    [citems{k}, sts1] = expand(citems{k}, eflag, tropts);
    sts = sts || sts1;
end;

item = subsasgn(item, ntgt_input, citems);
if sts && item.cfg_item.expanded == 0
    item.cfg_item.expanded = 1;
end;
