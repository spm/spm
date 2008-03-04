function [id stop val] = list(item, spec, tropts, fn)

% LIST function for cfg trees
% This function searches the cfg tree for certain entries.
%
% [id stop val] = list(item, spec, tropts[, fieldname])
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
% There are four pseudo-fieldnames which allow to obtain information useful
% to build e.g. a user interface for cfg trees:
% 'class' - returns the class of the current item
% 'level' - returns the level in the tree. Since data is collected
%           pre-order, children are listed after their parents. Identical
%           levels of subsequent nodes denote siblings, whereas decreasing
%           levels of subsequent nodes denote siblings of the parent node.
% 'all_set' - return all_set status of subtree rooted at item, regardless
%             whether list will descend into it or not
% 'all_set_item' - return all_set_item status of current node (i.e. whether
%                  all integrity conditions for this node are fulfilled)
%                  For in-tree nodes this can be different from all_set.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: list.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

if match(item, spec)
    id = {struct('type', {}, 'subs', {})};
    stop = false;
    if nargin > 3
        specialfn = {'class','level','all_set','all_set_item'};
        for k = 1:numel(fn)
            if any(strcmp(fn{k}, fieldnames(item)))
                val{k} = {subsref(item, substruct('.', fn{k}))};
            elseif any(strcmp(fn{k}, specialfn))
                switch fn{k}
                    case 'class'
                        val{k} = {class(item)};
                    case 'level'
                        val{k} = {tropts.clvl};
                    case 'all_set'
                        val{k} = {all_set(item)};
                    case 'all_set_item'
                        val{k} = {all_set_item(item)};
                end;
            else
                val{k} = {{}};
            end;
        end;
    else
        val = {};
    end;
else
    id = {};
    stop = [];
    if nargin > 3
        val = cell(size(fn));
        [val{:}] = deal({});        
    else
        val = {};
    end;
end;
