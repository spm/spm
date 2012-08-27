function [id, stop, val] = list(item, spec, tropts, fn)

% function [id, stop, val] = list(item, spec, tropts, fn)
% Find items in a cfg tree rooted at item that match a specification spec.
% By default, the filled configuration tree is searched (i.e. the
% val-branches of cfg_repeat and cfg_choice nodes). 
% See MATCH for help about spec data structure.
%
% Traversal options
% struct with fields
% stopspec - match spec to stop traversal
% dflag    - traverse val or values tree
% clvl     - current level in tree
% mlvl     - maximum level to traverse - range 1 (top level only) to
%            Inf (all levels)
% cnt      - #items found so far
% mcnt    - max #items to find
% List will stop descending into subtrees if one of the conditions
% following conditions are met: item matches stopspec, clvl >= mlvl, cnt >=
% mcnt. Flag stop is true for nodes where traversal has stopped
% (i.e. items where tropts has stopped further traversal).
%
% A cell list of subsref ids to matching nodes will be returned. The id of
% this node is returned before the id of its matching children.
% If the root node of the tree matches, the first id returned will be an
% empty substruct.
% If a cell list of fieldnames is given, then the contents of these fields
% will be returned in the cell array val. If one of the fields does not
% exist, a cell with an empty entry will be returned.
% There are five pseudo-fieldnames which allow to obtain information useful
% to build e.g. a user interface for cfg trees:
% 'class' - returns the class of the current item
% 'level' - returns the level in the tree. Since data is collected
%           pre-order, children are listed after their parents. Identical
%           levels of subsequent nodes denote siblings, whereas decreasing
%           levels of subsequent nodes denote siblings of the parent node.
% 'all_set' - return all_set status of subtree rooted at item, regardless
%             whether list will descend into it or not
% 'all_set_item' - return all_set_item status of current node (i.e. whether
%                  all integrity conditions for this node are fulfilled).
%                  For in-tree nodes this can be different from all_set.
% 'showdoc' - calls showmydoc to display the help text and option hints for
%             the current item (without recursive calls for .val/.values
%             items).
%
% This function is identical for all cfg_intree classes.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: list.m 4864 2012-08-27 13:57:31Z volkmar $

rev = '$Rev: 4864 $'; %#ok

if match(item, spec)
    id = {struct('type', {}, 'subs', {})};
    tropts.cnt = tropts.cnt + 1;
    stop = false;
    if nargin > 3
        val = cell(1,numel(fn));
        for k = 1:numel(fn)
            switch fn{k}
                case 'class'
                    val{k} = {class(item)};
                case 'level'
                    val{k} = {tropts.clvl};
                case 'all_set'
                    val{k} = {all_set(item)};
                case 'all_set_item'
                    val{k} = {all_set_item(item)};
                case 'showdoc'
                    val{k} = {showmydoc(item,'')};
                case fieldnames(item)
                    val{k} = {subsref(item, substruct('.', fn{k}))};
                otherwise
                    val{k} = {{}};
            end;
        end;
    else
        val = {};
    end;
    matched = true;
else
    id = {};
    stop = [];
    if nargin > 3
        val = cell(size(fn));
        [val{:}] = deal({});
    else
        val = {};
    end;
    matched = false;
end;

if (isempty(tropts.stopspec) || ~match(item, tropts.stopspec)) && (tropts.cnt < tropts.mcnt) && (tropts.clvl < tropts.mlvl)
    tname = treepart(item, tropts.dflag);
    % subscripts to child nodes
    idsubs = substruct('.', tname, '{}', {});
    citems = subsref(item, idsubs(1));
    tropts.clvl = tropts.clvl + 1;
    for k = 1:numel(citems)
        if nargin > 3
            [id1 stop1 val1] = list(citems{k}, spec, tropts, fn);
            for l = 1:numel(fn)
                val{l} = [val{l}(:); val1{l}(:)]';
            end;
        else
            [id1 stop1] = list(citems{k}, spec, tropts);
        end;
        if ~isempty(id1)
            idsubs(2).subs = {k};
            for l = 1:numel(id1)
                id = [id(:); {[idsubs id1{l}]}]';
            end;
            stop = [stop stop1];
            tropts.cnt = tropts.cnt + numel(id1);
            if tropts.cnt >= tropts.mcnt
                if matched
                    stop(1) = true;
                end;
                break;
            end;
        end;
    end;
else
    if matched
        stop = true;
    end;
end;
